import sys,pathlib
import numpy as np
import h5py as h5
import uproot as ur
from typing import Dict, Any, Tuple, Optional
from abc import ABC, abstractmethod
# Needed for the test function
# import os
# this_dir = os.path.dirname(os.path.abspath(__file__))
# sys.path.append(str(pathlib.Path('{}/../'.format(this_dir)).absolute()))
from util.calcs import embed_array

# Classes for data buffers, to be used by the post-processors (such as JetFinder).
# NOTE: Our HepMC -> HDF5 conversion in conversion.py also uses some buffering logic,
#       but it doesn't implement this class. Maybe it can eventually be updated? - Jan

class BufferFlushHandler(ABC):
    """Abstract base class for handling buffer flush operations."""

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

class DummyFlushHandler(BufferFlushHandler):
    """Example flush handler that prints what would be written to file."""

    def __init__(self, filename: str):
        self.filename = filename
        self.verbose = True

    def flush(self, data: Dict[str, np.ndarray], start_event: int, end_event: int, nevents: int):
        if(self.verbose):
            print("\tFlushing events {}-{} to {}".format(start_event,end_event-1,self.filename))
            for key, array in data.items():
                print("  {}: shape {}, dtype {}".format(key,array.shape,array.dtype))

class HDF5FlushHandler(BufferFlushHandler):
    """Flushes data to an HDF5 file.."""

    def __init__(self, filename: str):
        self.filename = filename
        self.status = 'w'
        self.copts = 9
        self.verbose = False
        self.f = None

    def _set_status(self):
        """
        Automatically determines if we need to write the file,
        or if it already exists and we're appending.
        """
        if(pathlib.Path(self.filename).exists()):
            self.status = 'a'

    def flush(self, data: Dict[str, np.ndarray], start_event: int, end_event: int, nevents: int):
        self._set_status()
        if(self.verbose): print("\tFlushing events {}-{} to {}".format(start_event,end_event-1,self.filename))
        self.f = h5.File(self.filename,self.status)
        for key, array in data.items():

            if((self.status == 'w') or (key not in self.f.keys())):
                if(self.verbose): print('\tCreating dset {}'.format(key))
                dset = self._create_dataset(key,array,nevents)
            else:
                if(self.verbose): print('\tLoading dset {}'.format(key))
                dset = self.f[key]

            dset[start_event:end_event] = array
        self.f.close()
        self.f = None

    def _create_dataset(self,key,array,nevents):
        # Need to initialize the dataset in the HDF5 file, which requires
        # getting the right shape (namely the 1st dimension!).
        dset_shape = tuple((nevents,) + array.shape[1:])
        dset = self.f.create_dataset(key, data=np.zeros(dset_shape,dtype=array.dtype),compression='gzip',compression_opts=self.copts)
        return dset

class BufferArray:
    """
    A wrapper around numpy arrays that handles circular buffer indexing automatically.
    """

    def __init__(self, array: np.ndarray, buffer: 'Buffer'):
        self._array = array
        self._buffer = buffer

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
            # Single event index - handle circular buffer logic
            buffer_position = key % self._buffer.buffer_size
            self._array[buffer_position] = value
            self._buffer._written[self._buffer._key][buffer_position] = True
            self._buffer._number_written[self._buffer._key] = np.sum(self._buffer._written[self._buffer._key])
            self._buffer._total_events_processed = max(self._buffer._total_events_processed, key + 1)

            # Track that this event position has been written to
            # NOTE: This buffer code is currently too contrived
            self._buffer._current_event_positions.add(buffer_position)
            self._buffer._current_size = len(self._buffer._current_event_positions)
            self._buffer._check_and_flush_if_needed(key)

        elif isinstance(key, tuple):
            # Multi-dimensional indexing like [event_index, i, :j]
            first_idx = key[0]
            rest_idx = key[1:]

            if isinstance(first_idx, int):

                # First dimension is an event index - apply circular buffer logic
                buffer_position = first_idx % self._buffer.buffer_size
                self._array[(buffer_position,) + rest_idx] = value
                self._buffer._written[self._buffer._key][buffer_position] = True
                self._buffer._number_written[self._buffer._key] = np.sum(self._buffer._written[self._buffer._key])
                self._buffer._total_events_processed = max(self._buffer._total_events_processed, first_idx + 1)

                # Track that this event position has been written to
                # NOTE: This buffer code is currently too contrived
                self._buffer._current_event_positions.add(buffer_position)
                self._buffer._current_size = len(self._buffer._current_event_positions)
                self._buffer._check_and_flush_if_needed(first_idx)

            else:
                raise ValueError("Slicing not understood or implemented in this way.")
                # First dimension is a slice/fancy index - pass through directly
                self._array[key] = value
                self._buffer._number_written[self._buffer._key] += len(key[0]) # no idea if this works -- will probably never call it

        else:
            raise ValueError("Slicing not understood or implemented in this way.")
            # Single slice or other indexing on first dimension - pass through directly
            self._array[key] = value
            self._buffer._number_written[self._buffer._key] += len(key[0]) # no idea if this works -- will probably never call it

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

class Buffer:
    """
    A dictionary-like buffer that maintains fixed-size numpy arrays and
    automatically flushes data when the buffer fills up.
    """

    def __init__(self, buffer_size: int = 100, filename=None, flush_handler: Optional[BufferFlushHandler] = None):
        """
        Initialize the buffer.

        Args:
            buffer_size: Maximum number of events to store before flushing
            flush_handler: Handler for flush operations (optional)
        """
        self.filename = filename
        self.buffer_size = buffer_size
        self.nevents = -1
        self.flush_handler = flush_handler
        if(self.flush_handler is None):
            self.flush_handler = HDF5FlushHandler(self.filename)

        # Internal storage
        self._buffer_arrays: Dict[str, BufferArray] = {}
        self._array_specs: Dict[str, Tuple] = {}  # Store (shape, dtype) for each array

        # Track buffer state - much simpler approach
        self._current_size = 0  # Number of events currently in buffer
        self._total_events_processed = 0
        self._buffer_start_event = 0
        self._current_event_positions = set()  # Track which event positions have been written to

        self._written : Dict[str, np.ndarray] = {}
        self._number_written: Dict[str, int] = {}
        self._key = None

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

        # if(self.nevents > -1):
        #     self.nevents = shape[0]
        # else:
        #     print('self.nevents = {}, shape[0] = {}'.format(self.nevents,shape[0]))
        #     assert self.nevents == shape[0]

    def _check_and_flush_if_needed(self, event_index: int):
        """Check if we need to flush before processing this event."""
        # buffer_position = event_index % self.buffer_size

        do_flush = True
        # print('Check flush')
        for key,val in self._number_written.items():
            # print('\t-> {}, {}'.format(key,val))
            if(val != self.buffer_size):
                do_flush = False
                # print('\t\t->False')
                break

        if(do_flush):
            self._flush_buffer()

    def _flush_buffer(self):
        """Flush the current buffer contents."""
        if self.flush_handler and len(self._current_event_positions) > 0:
            current_size = len(self._current_event_positions)
            # Create a view of only the filled portion of each array
            flush_data = {}
            for key, buffer_array in self._buffer_arrays.items():
                flush_data[key] = buffer_array._array[:current_size].copy()

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
        """Set an entire array in the buffer (not typically used for event-by-event)."""
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

    def set(self,key: str, index: int, value: Any):
        """
        For setting values of entries in the Buffer. (i.e. particular BufferArrays).
        Handles embedding/zero-padding as needed.
        In general, this is the function one should use for putting data into the buffer.
        """
        if(isinstance(value,int) or isinstance(value,float)):
            self[key][index] = value
        else:
            if(not isinstance(value,np.ndarray)):
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
            shape: Shape including event dimension (e.g., (max_events, n_jets, 4))
            dtype: Data type for the array
        """
        if(isinstance(shape,int)):
            shape = (shape,)
        if key not in self._buffer_arrays:
            arr_shape = (1,) + shape
            self._initialize_array(key, arr_shape, dtype)
        return self._buffer_arrays[key]

# Example usage
def main(args):

    # Create a buffer with file flush handler
    flush_handler = HDF5FlushHandler("buffer_test.h5")
    buffer = Buffer(buffer_size=3, flush_handler=flush_handler)

    # Create arrays in buffer with some known maximum dimensionss.
    buffer.create_array('jets.N', dtype=np.int32)  # Scalar per event
    buffer.create_array('jets.Pmu', (10, 4), np.float64)  # Up to 10 jets, 4-momentum each
    buffer.create_array('event.weight', dtype=np.float64)  # Event weight

    # Process events
    nevents = 8
    buffer.SetNEvents(nevents)
    for event_index in range(nevents):
        # Simulate jet data for this event
        n_jets = np.random.randint(2, 6)
        jet_vectors = np.random.randn(n_jets, 4)  # (N_jets, 4-momentum)
        jet_ordering = np.arange(n_jets)  # Dummy ordering
        event_weight = np.random.uniform(0.5, 2.0)

        print("Processing event {} (buffer position {})".format(event_index, event_index % buffer.buffer_size))

        # Fill all branches for this event

        buffer.set('jets.N',event_index,len(jet_vectors))
        buffer.set('jets.Pmu',event_index,np.vstack([jet_vectors[i] for i in jet_ordering]))
        buffer.set('event.weight',event_index,event_weight)

        # buffer['jets.N'][event_index] = len(jet_vectors)

        # # Filling this is a bit verbose -- would be nice to make a function for this,
        # # to hide some of the details. - Jan
        # buffer['jets.Pmu'][event_index] = embed_array(np.vstack([jet_vectors[i] for i in jet_ordering]),buffer['jets.Pmu'][event_index].shape)

        # # embed_array_inplace(np.vstack([jet_vectors[i] for i in jet_ordering]),
        # #                    buffer['jets.Pmu'][event_index])

        # buffer['event.weight'][event_index] = event_weight

        print("  -> Jets: {}, buffer current size: {}".format(n_jets,buffer._current_size))
        print()

    # Final flush for remaining data
    print("Final flush of remaining {} events:".format(buffer._current_size))
    buffer.flush()
    print("\nBuffer info:", buffer.get_buffer_info())

if __name__ == "__main__":
    main(sys.argv)


########################################################
# Uproot/dask-related stuff, for use with Delphes reading.
########################################################

class IndexableLazyLoader:
    def __init__(self, files, tree_path, expressions):
        self.files = files
        self.tree_path = tree_path
        self.expressions = expressions
        self._cached_branches = {}
        self._total_entries = None

        # Get field info by checking what actually exists
        with ur.open(files[0]) as f:
            available_expressions = [expr for expr in expressions if expr in f[tree_path]]
        self.fields = available_expressions

    def __getitem__(self, branch_name):
        """Return a branch that can be sliced"""
        if branch_name not in self._cached_branches:
            self._cached_branches[branch_name] = LazyBranch(self.files, self.tree_path, branch_name)
        return self._cached_branches[branch_name]

    def __len__(self):
        """Return total number of entries across all files"""
        if self._total_entries is None:
            total = 0
            for file_path in self.files:
                with ur.open(file_path) as f:
                    total += f[self.tree_path].num_entries
            self._total_entries = total
        return self._total_entries

class LazyBranch:
    def __init__(self, files, tree_path, branch_name):
        self.files = files
        self.tree_path = tree_path
        self.branch_name = branch_name
        self._cached_data = None

    def __getitem__(self, slice_obj):
        """Handle slicing like delphes_arr[branch][start:stop]"""
        if self._cached_data is None:
            # Load the entire branch when first accessed
            self._cached_data = ur.concatenate(self.files,
                                             expressions=[self.branch_name],
                                             tree_path=self.tree_path,
                                             library="ak")[self.branch_name]
        return self._cached_data[slice_obj]

    def __len__(self):
        """Return length of this branch (same as total entries)"""
        if self._cached_data is not None:
            return len(self._cached_data)
        else:
            # Calculate without loading the data
            total = 0
            for file_path in self.files:
                with ur.open(file_path) as f:
                    total += f[self.tree_path].num_entries
            return total