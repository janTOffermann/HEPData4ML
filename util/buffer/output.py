import pathlib
import numpy as np
import h5py as h5
import ROOT as rt
import subprocess as sub
from typing import Dict, Any, Tuple, Optional, Union
from abc import ABC, abstractmethod
from util.math.embedding import embed_array
from util.hdf5.hdf5 import MergeH5
from util.misc.timing import profile_method, profile_block

# Classes for output data buffers, to be used by the post-processors (such as JetFinder).
# NOTE: Our HepMC -> HDF5 conversion in conversion.py also uses some buffering logic,
#       but it doesn't implement this class. Maybe it can eventually be updated? - Jan

class BufferFlushHandler(ABC):
    """Abstract base class for handling buffer flush operations."""

    def __init__(self):
        self.print_prefix = 'BufferFlushHandler:'

    @abstractmethod
    def flush(self, data: Dict[str, np.ndarray], start_event: int, end_event: int):
        """
        Handle flushing of buffer data.

        Args:
            data: Dictionary of numpy arrays to flush
            start_event: Starting event index (inclusive)
            end_event: Ending event index (exclusive)
        """
        pass

    def SetFilename(self, filename: str):
        self.filename = filename

    def SetVerbosity(self, verbose: bool):
        self.verbose = verbose

    def _print(self,val:str):
        print('{} {}'.format(self.print_prefix,val))

class DummyFlushHandler(BufferFlushHandler):
    """Example flush handler that prints what would be written to file."""

    def __init__(self, filename: str):
        self.filename = filename
        self.verbose = True
        self.print_prefix = 'DummyFlushHandler:'

    def flush(self, data: Dict[str, np.ndarray], start_event: int, end_event: int, nevents: int):
        if(self.verbose):
            self._print("\tFlushing events {}-{} to {}".format(start_event,end_event-1,self.filename))
            for key, array in data.items():
                self._print("  {}: shape {}, dtype {}".format(key,array.shape,array.dtype))

class HDF5FlushHandler(BufferFlushHandler):
    """Flushes data to an HDF5 file.."""

    def __init__(self, filename: str, copts=0):
        self.filename = filename
        self.status = 'w'
        self.copts = copts
        self.verbose = False
        self.f = None
        self.buffer_size = None # for chunking storage
        self.print_prefix = 'HDF5FlushHandler:'

    def _init_file(self): # Note: Keeping the file open was an attempt to speed things up, not clear that it makes a difference. - Jan
        self._set_status()
        self.f = h5.File(self.filename,self.status)

    def SetCompression(self,val:int):
        if(val < 0): val = 0
        elif(val > 9): val = 9
        self.copts = val

    def SetBufferSize(self,val:int):
        self.buffer_size = val

    def _set_status(self):
        """
        Automatically determines if we need to write the file,
        or if it already exists and we're appending.
        """
        if(pathlib.Path(self.filename).exists()):
            self.status = 'a'

    def flush(self, data: Dict[str, np.ndarray], start_event: int, end_event: int, nevents: int):

        if(self.f is None):
            self._init_file()

        if(self.verbose): self._print("\tFlushing events {}-{} to {}".format(start_event,end_event-1,self.filename))
        # self.f = h5.File(self.filename,self.status)
        for key, array in data.items():

            if(key not in self.f.keys()):
                if(self.verbose): self._print('\tCreating dset {}'.format(key))
                dset = self._create_dataset(key,array,nevents)
            else:
                if(self.verbose): self._print('\tLoading dset {}'.format(key))
                dset = self.f[key]

            dset[start_event:end_event] = array
        # self.f.close()
        # self.f = None

    def _create_dataset(self,key,array,nevents):
        # Need to initialize the dataset in the HDF5 file, which requires
        # getting the right shape (namely the 1st dimension!).
        dset_shape = tuple((nevents,) + array.shape[1:])
        if(self.buffer_size is not None):
            chunks = tuple((self.buffer_size,) + array.shape[1:])
        else:
            chunks = True # or set to False?
        dset = self.f.create_dataset(
            key,
            data=np.zeros(dset_shape,dtype=array.dtype),
            compression='gzip',
            compression_opts=self.copts,
            chunks=chunks
        )
        return dset

    def close(self,output_file:Optional[str]=None, copts=0):

        # First, close the buffer output.
        self.f.close()

        # Now, we optionally merge the buffer output
        # into the  provided output_file.
        if(output_file is not None):
            MergeH5(output_file,
                    self.filename,
                    copts=copts,
                    delete_input=True
            )
        return

class BufferArray:
    """
    A wrapper around numpy arrays that handles circular buffer indexing automatically.
    """

    def __init__(self, array: np.ndarray, buffer: 'OutputBuffer'):
        self._array = array
        self._buffer = buffer
        self.print_prefix = 'BufferArray:'

    def __getitem__(self, key):
        """Get data, automatically handling circular buffer indexing for the first dimension."""
        if isinstance(key, int):
            # Single event index - handle circular buffer logic
            buffer_position = key % self._buffer.buffer_size
            return self._array[buffer_position]

        elif isinstance(key, tuple):
            # Multi-dimensional indexing like [event_index, i, :j]
            first_idx = key[0]
            rest_idx = key[1:]

            if isinstance(first_idx, int):
                # First dimension is an event index - apply circular buffer logic
                buffer_position = first_idx % self._buffer.buffer_size
                return self._array[(buffer_position,) + rest_idx]
            else:
                # First dimension is a slice/fancy index - pass through directly
                # This handles cases like buffer['jets.Pmu'][:10, i, j]
                return self._array[key]

        else:
            # Single slice or other indexing on first dimension - pass through directly
            # This handles cases like buffer['jets.Pmu'][:10] or buffer['jets.Pmu'][::2]
            return self._array[key]

    def __setitem__(self, key, value):
        """Set data, automatically handling circular buffer indexing for the first dimension."""
        if isinstance(key, int):

            # TODO: This code needs some cleaning up; it has been modified a couple times, and
            #       is probably unnecessarily convoluted.
            buffer_position = key % self._buffer.buffer_size
            self._buffer._written[self._buffer._key][buffer_position] = True
            self._buffer._number_written[self._buffer._key] = np.sum(self._buffer._written[self._buffer._key])
            self._buffer._total_events_processed = max(self._buffer._total_events_processed, key + 1)

            # Track that this event position has been written to
            self._buffer._current_event_positions.add(buffer_position)
            self._buffer._current_size = len(self._buffer._current_event_positions)
            self._buffer._check_and_flush_if_needed(key)
            self._array[buffer_position] = value # after the possible flush, it's safe to actually modify self._array


        elif isinstance(key, tuple):
            # Multi-dimensional indexing like [event_index, i, :j]
            first_idx = key[0]
            rest_idx = key[1:]

            if isinstance(first_idx, int):

                buffer_position = first_idx % self._buffer.buffer_size
                self._buffer._written[self._buffer._key][buffer_position] = True
                self._buffer._number_written[self._buffer._key] = np.sum(self._buffer._written[self._buffer._key])
                self._buffer._total_events_processed = max(self._buffer._total_events_processed, first_idx + 1)

                # Track that this event position has been written to
                self._buffer._current_event_positions.add(buffer_position)
                self._buffer._current_size = len(self._buffer._current_event_positions)
                self._buffer._check_and_flush_if_needed(first_idx)
                self._array[(buffer_position,) + rest_idx] = value # after the possible flush, it's safe to actually modify self._array

            else:
                raise ValueError("Slicing not understood or implemented in this way.")
                # First dimension is a slice/fancy index - pass through directly
                # self._array[key] = value
                # self._buffer._number_written[self._buffer._key] += len(key[0]) # no idea if this works -- will probably never call it

        else:
            raise ValueError("Slicing not understood or implemented in this way.")
            # Single slice or other indexing on first dimension - pass through directly
            # self._array[key] = value
            # self._buffer._number_written[self._buffer._key] += len(key[0]) # no idea if this works -- will probably never call it

    @property
    def shape(self):
        """Return the shape of the underlying array."""
        return self._array.shape

    @property
    def dtype(self):
        """Return the dtype of the underlying array."""
        return self._array.dtype

    def __repr__(self):
        return f"BufferArray(shape={self.shape}, dtype={self.dtype})"

    def _print(self,val:str):
        print('{} {}'.format(self.print_prefix,val))

class OutputBuffer:
    """
    A dictionary-like buffer that maintains fixed-size numpy arrays and
    automatically flushes data when the buffer fills up.
    """

    def __init__(self, buffer_size: int = 100, filename:Optional[str]=None, flush_handler: Optional[BufferFlushHandler] = None, copts=9):
        """
        Initialize the buffer.

        Args:
            buffer_size: Maximum number of events to store before flushing
            flush_handler: Handler for flush operations (optional)
        """
        self.filename = filename
        self.nevents = -1
        self.flush_handler = flush_handler
        if(self.flush_handler is None):
            self.flush_handler = HDF5FlushHandler(self.filename,copts)

        # Internal storage
        self._buffer_arrays: Dict[str, BufferArray] = {}
        self._array_specs: Dict[str, Tuple] = {}  # Store (shape, dtype) for each array

        # Track buffer state
        self._current_size = 0  # Number of events currently in buffer
        self._total_events_processed = 0
        self._buffer_start_event = 0
        self._current_event_positions = set()  # Track which event positions have been written to

        self._written : Dict[str, np.ndarray] = {}
        self._number_written: Dict[str, int] = {}
        self._key = None

        self.SetBufferSize(buffer_size)

        self.print_prefix = 'OutputBuffer:'

    def SetBufferSize(self,buffer_size:int):
        self.buffer_size = buffer_size
        if(self.flush_handler is not None):
            self.flush_handler.SetBufferSize(self.buffer_size)
        return

    def SetFilename(self, filename: str):
        self.filename = filename
        self.flush_handler.SetFilename(self.filename)

    def SetNEvents(self,nevents: int):
        self.nevents = nevents

    def _initialize_array(self, key: str, shape: Tuple, dtype: np.dtype):
        """Initialize a new array in the buffer with the given specifications."""
        full_shape = (self.buffer_size,) + shape[1:]  # Replace first dim with buffer_size
        array = np.zeros(full_shape, dtype=dtype)
        self._buffer_arrays[key] = BufferArray(array, self)
        self._array_specs[key] = (shape, dtype)
        self._written[key] = np.full(self.buffer_size,False)
        self._number_written[key] = 0

    @profile_method('OutputBuffer._check_and_flush_if_needed')
    def _check_and_flush_if_needed(self, event_index: int):
        """Check if we need to flush before processing this event."""
        buffer_position = event_index % self.buffer_size
        if(not ((buffer_position == 0) and (event_index != 0))):
            return

        for key,val in self._number_written.items():
            if(val != self.buffer_size):
                return

        self._flush_buffer()

    @profile_method('OutputBuffer._flush_buffer')
    def _flush_buffer(self):
        """Flush the current buffer contents."""
        if self.flush_handler and len(self._current_event_positions) > 0:
            current_size = len(self._current_event_positions)
            # Create a view of only the filled portion of each array
            flush_data = {}
            for key, buffer_array in self._buffer_arrays.items():

                flush_data[key] = buffer_array._array[:current_size].copy() # removed the ndarray.copy() function

            self.flush_handler.flush(
                flush_data,
                self._buffer_start_event,
                self._buffer_start_event + current_size,
                self.nevents
            )

        # Reset buffer state
        self._buffer_start_event += len(self._current_event_positions)
        self._current_event_positions.clear()
        self._current_size = 0
        self._written = {key:np.full(self.buffer_size,False) for key in self._written.keys()}
        self._number_written = {key:0 for key in self._number_written.keys()}

    @property
    def _current_size_prop(self):
        """Current number of filled positions in the buffer (for compatibility)."""
        return self._current_size

    def __getitem__(self, key: str) -> BufferArray:
        """Get a BufferArray from the buffer (dictionary-like access)."""

        if key not in self._buffer_arrays:
            raise KeyError(f"Key '{key}' not found in buffer")
        self._key = key
        return self._buffer_arrays[key]

    def __setitem__(self, key: str, value: np.ndarray):
        """Set an entire array in the buffer (not used for event loop approach)."""
        if not isinstance(value, np.ndarray):
            value = np.array(value)
        self._key = key

        if key not in self._buffer_arrays:
            self._initialize_array(key, value.shape, value.dtype)

        # Ensure the array fits in the buffer
        copy_size = min(value.shape[0], self.buffer_size)
        self._buffer_arrays[key]._array[:copy_size] = value[:copy_size]
        self._current_size = max(self._current_size, copy_size)
        for i in range(copy_size):
            self._written[key][i] = True
        self._number_written[key] = np.sum(self._written[key])

    @profile_method('OutputBuffer.set')
    def set(self,key: str, index: Optional[int], value: Any):
        """
        For setting values of entries in the Buffer. (i.e. particular BufferArrays).
        Handles embedding/zero-padding as needed.
        In general, this is the function one should use for putting data into the buffer.
        """

        index = slice(None) if index is None else index

        if(isinstance(value,int) or isinstance(value,float)):
            self[key][index] = value
        else:
            if(isinstance(value,list)):
                value = np.array(value)
            self[key][index] = embed_array(value,self[key][index].shape)
        return

    def __contains__(self, key: str) -> bool:
        """Check if a key exists in the buffer."""
        return key in self._buffer_arrays

    def keys(self):
        """Return the keys in the buffer."""
        return self._buffer_arrays.keys()

    def flush(self):
        """Manually flush the current buffer contents."""
        self._flush_buffer()

    def get_buffer_info(self) -> Dict[str, Any]:
        """Get information about the current buffer state."""
        return {
            'buffer_size': self.buffer_size,
            'current_size': self._current_size,
            'total_events_processed': self._total_events_processed,
            'buffer_start_event': self._buffer_start_event,
            'arrays': {key: {'shape': arr.shape, 'dtype': arr.dtype}
                      for key, arr in self._buffer_arrays.items()}
        }

    def create_array(self, key: str, shape: Tuple=(), dtype: np.dtype = np.float64):
        """
        Explicitly create an array in the buffer with specified shape and dtype.

        Args:
            key: Array key
            shape: Shape excluding event dimension
            dtype: Data type for the array
        """
        if(isinstance(shape,int)):
            shape = (shape,)
        if key not in self._buffer_arrays:
            arr_shape = (1,) + shape
            self._initialize_array(key, arr_shape, dtype)
        return self._buffer_arrays[key]

    def close(self,output_file:Optional[str]=None):
        self.flush_handler.close(output_file)

    def _print(self,val:str):
        print('{} {}'.format(self.print_prefix,val))


class RootOutputBuffer:
    """
    Similar in usage to OutputBuffer, but doesn't do complex buffering;
    it writes output to a ROOT file in the "standard" way,
    """

    def __init__(self, filename:Optional[str]=None, tree_name:Optional[str]=None):
        self.filename = None
        self.tree_name = tree_name if tree_name is not None else "hepdata4ml_tree"
        self.print_prefix = 'RootOutputBuffer:'
        self.buffers = {} # each buffer will be of length 1 w.r.t. number of events
        self.buffer_is_vector = {}
        self.n_filled = 0 # keep track of what buffers have been filled
        self.is_scalar = {}
        self.f = None
        self.t = None
        self.clone_tree = None
        self.init_status = False

        self.SetFilename(filename)

    def SetFilename(self, filename: str):
        if(filename is None):
            return
        self.filename = filename

    def SetTreeName(self,tree_name:str):
        self.tree_name = tree_name

    def SetCloneTree(self,tree:rt.TTree):
        self.clone_tree = tree

    def _init_tree(self,in_tree=None):
        self.f = rt.TFile(self.filename,'RECREATE')

        if(in_tree is None and self.clone_tree is None):
            self.t = rt.TTree(self.tree_name,self.tree_name)
        elif(self.clone_tree is not None):
            self.t = self.clone_tree.CloneTree(0) # clones the TTree structure, does *not* copy over any entries
        else:
            self.t = in_tree.CloneTree(0) # clones the TTree structure, does *not* copy over any entries
        self.init_status = True
        return

    def GetTreeName(self):
        return self.tree_name

    def create_array(self, key: str, shape: Tuple=(), dtype: np.dtype = np.float64):
        """
        Explicitly create an array in the buffer with specified shape and dtype.
        Based on OutputBuffer.create_array()
        Args:
            key: Array key
            shape: Shape excluding event dimension (used to determine if scalar or std::vector)
            dtype: Data type for the array
        """
        if(key in self.buffers.keys()):
            return

        if(not self.init_status):
            self._init_tree()

        if(shape==()): # scalar -- one per event
            self._init_scalar_branch(key,dtype)
        else:
            self._init_vector_branch(key,shape, dtype)
        # self.n_filled[key] = 0
        return

    def _init_scalar_branch(self,key,dtype):
        self.buffers[key] = np.zeros(1,dtype=dtype)
        self.buffer_is_vector[key] = False

        #TODO: Support more types?
        if(dtype == np.dtype('float')):
            self.t.Branch(key,self.buffers[key],'{}/D'.format(key))
        elif(dtype == np.dtype('i4')):
            self.t.Branch(key,self.buffers[key],'{}/I'.format(key))
        elif(dtype == np.dtype('uint')):
            self.t.Branch(key,self.buffers[key],'{}/i'.format(key))
        elif(dtype == np.dtype('i8')): # also covers long
            self.t.Branch(key,self.buffers[key],'{}/L'.format(key))
        elif(dtype == np.dtype('ulong')): # also covers long
            self.t.Branch(key,self.buffers[key],'{}/l'.format(key))
        elif(dtype == np.dtype('short')):
            self.t.Branch(key,self.buffers[key],'{}/S'.format(key))
        elif(dtype == np.dtype('ushort')):
            self.t.Branch(key,self.buffers[key],'{}/s'.format(key))
        elif(dtype == np.dtype('bool')):
            self.t.Branch(key,self.buffers[key],'{}/o'.format(key))
        else: # not recognized
            self._print('Warning: dtype {} not recognized for branch {}.'.format(dtype,key))
        return

    def _init_vector_branch(self,key,shape, dtype):
        dtype_str = 'double'
        for type_str in ['int','uint','short','ushort','bool']:
            if(dtype == np.dtype(type_str)):
                dtype_str = type_str
                break
        if(dtype == np.dtype('int32')): # special case not caught above TODO: Running into some issues when pushing back to vector<vector<int>> in this case?
            dtype_str = 'int'

        # For now, we will support 1D, 2D and 3D vectors
        if(len(shape) == 1):
            self.buffers[key] = rt.std.vector[dtype_str]()
        elif(len(shape) == 2):
            self.buffers[key] = rt.std.vector[rt.std.vector[dtype_str]]()
        elif(len(shape) == 3):
            self.buffers[key] = rt.std.vector[rt.std.vector[rt.std.vector[dtype_str]]]()
        else:
            self._print('Warning: vector branch of dimension {} not supported.'.format(len(shape)))
            return
        self.buffer_is_vector[key] = True

        self.t.Branch(key,self.buffers[key])
        return

    def set(self,key: str, index: Optional[int], value: Any):
        """
        For setting values of entries in the buffer.
        In general, this is the function one should use for putting data into the buffer.
        """
        # We need to cover multiple cases: scalar-type branches, and vector-type branches.
        # For the case of vectors, they can be multi-dimensional (e.g. vector<vector<Double_t>>).
        try:
            is_scalar = self.is_scalar[key]
        except:
            # need to determine if this is a scalar branch or not
            is_scalar = False
            if(isinstance(value,int) or isinstance(value,float)):
               is_scalar = True
            self.is_scalar[key] = is_scalar

        # Check if we need to flush the buffer. We do this if we find that
        # we're filling the next event for this buffer, as identified by
        # the `index` argument ; this works as long
        # as the code that's leveraging this class is filling all the buffers
        # for a single event before moving on to the next one.
        # (which is a pretty sensible assumption) - Jan
        event_index = index
        if(isinstance(index,tuple)):
            event_index = index[0]

        if(self.n_filled < event_index):
            self.flush()

        if(is_scalar):
            self.buffers[key][0] = value # buffer is a 1D length-1 array
        else: # non-scalar -- this possibly gets more complex

            if(isinstance(index,tuple)): # slicing: gets a bit complex

                if(len(index) > 2):
                    self._print('Warning: Indexing beyond 2 dims not (yet) supported for RootOutputBuffer.set().')
                    return

                if isinstance(value, np.ndarray) and not value.flags['C_CONTIGUOUS']:
                    value = np.ascontiguousarray(value)

                # We need to determine if we're filling the "next" entry in this vector,
                # or are overwriting an existing entry (unlikely!) or writing non-sequentially.
                # I should emphasize that these latter cases are highly unlikely given how
                # this class is expected to be used. - Jan
                current_length = self.buffers[key].size()
                if(index[1] == current_length):
                    self.buffers[key].push_back(value)

                elif(index[1] > current_length):
                    for i in range(index[1] - current_length):
                        self.buffers[key].push_back({}) # need to put some empty entries to pre-pad
                else:
                    self.buffers[key][index[1]] = value # reassign

            else:
                self.buffers[key].assign(value) # TODO: Does this work as expected?
        return

    def flush(self):
        self.t.Fill()
        self.n_filled += 1
        self._clear_buffers()

    def _clear_buffers(self):
        for key in self.buffers.keys():
            if(self.buffer_is_vector[key]):
                self.buffers[key].clear()


    def close(self,output_file:Optional[str]=None):
        self.t.Write()
        self.f.Close()

        if(output_file is not None):
            command = ['mv',self.filename,output_file]
            sub.check_call(command)
        return

    def keys(self):
        """Return the keys in the buffer."""
        return self.buffers.keys()

    def _print(self,val:str):
        print('{} {}'.format(self.print_prefix,val))