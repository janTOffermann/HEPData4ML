import pathlib
import numpy as np
import h5py as h5
from typing import Dict, Any, Tuple, Optional
from abc import ABC, abstractmethod

def embed_array_inplace(array, target, padding_value=0):
    """
    A generic function for embedding an input array into some target array.
    The array will be truncated or padded as needed.
    Modifies target in-place.
    """
    array_np = np.array(array)

    # Calculate the effective shape (minimum dimensions)
    effective_shape = tuple(min(s, t) for s, t in zip(array_np.shape, target.shape))

    # Create slices for both arrays
    slices = tuple(slice(0, dim) for dim in effective_shape)

    # Reset target to padding value
    target.fill(padding_value)

    # Copy data from array to target
    target[slices] = array_np[slices]
    return

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

class DummyFlushHandler(BufferFlushHandler):
    """Example flush handler that prints what would be written to file."""

    def __init__(self, filename: str):
        self.filename = filename

    def flush(self, data: Dict[str, np.ndarray], start_event: int, end_event: int, nevents: int):
        print(f"Flushing events {start_event}-{end_event-1} to {self.filename}")
        for key, array in data.items():
            print(f"  {key}: shape {array.shape}, dtype {array.dtype}")
        # Here you would implement actual file writing (HDF5, ROOT, etc.)

class HDF5FlushHandler(BufferFlushHandler):
    """Flushes data to an HDF5 file.."""

    def __init__(self, filename: str):
        self.filename = filename
        self.status = 'w'
        self.copts = 9

    def _set_status(self):
        if(pathlib.Path(self.filename).exists()):
            self.status = 'a'

    def flush(self, data: Dict[str, np.ndarray], start_event: int, end_event: int, nevents: int):
        self._set_status()
        print(f"Flushing events {start_event}-{end_event-1} to {self.filename}")
        f = h5.File(self.filename,self.status)
        for key, array in data.items():

            if(self.status == 'w'):
                # Need to initialize the dataset in the HDF5 file, which requires
                # getting the right shape (namely the 1st dimension!).
                dset_shape = (nevents) + array.shape[1:]
                dset = f.create_dataset(key, np.zeros(dset_shape), array.dtype,compression='gzip',compression_opts=self.copts)
            else:
                dset = f[key]
                dset[start_event:end_event] = array

            print(f"  {key}: shape {array.shape}, dtype {array.dtype}")
        f.close()
        # Here you would implement actual file writing (HDF5, ROOT, etc.)


class BufferArray:
    """
    A wrapper around numpy arrays that handles circular buffer indexing automatically.
    """

    def __init__(self, array: np.ndarray, buffer: 'Buffer'):
        self._array = array
        self._buffer = buffer

    def __getitem__(self, key):
        """Get data, automatically handling circular buffer indexing for the first dimension."""
        # NOTE: A little hacky, but we adjust self._buffer._number_written, even here in __getitem__.
        #       This is to deal with cases where values are adjusted in-place, with the side-effect that
        #       there will be unwanted behavior if one ever accesses BufferArray without writing to it.
        if isinstance(key, int):
            # Single event index - handle circular buffer logic
            # self._buffer._check_and_flush_if_needed(key)
            buffer_position = key % self._buffer.buffer_size
            self._buffer._number_written[self._buffer._key] += 1
            self._buffer._check_and_flush_if_needed(key)
            return self._array[buffer_position]

        elif isinstance(key, tuple):
            # Multi-dimensional indexing like [event_index, i, :j]
            first_idx = key[0]
            rest_idx = key[1:]

            if isinstance(first_idx, int):
                # First dimension is an event index - apply circular buffer logic
                self._buffer._check_and_flush_if_needed(first_idx)
                buffer_position = first_idx % self._buffer.buffer_size
                self._buffer._number_written[self._buffer._key] += 1
                return self._array[(buffer_position,) + rest_idx]
            else:
                # First dimension is a slice/fancy index - pass through directly
                # This handles cases like buffer['jets.Pmu'][:10, i, j]
                self._buffer._number_written[self._buffer._key] += len(key[0]) # TODO: Not sure if this works

                return self._array[key]

        else:
            # Single slice or other indexing on first dimension - pass through directly
            # This handles cases like buffer['jets.Pmu'][:10] or buffer['jets.Pmu'][::2]
            return self._array[key]

    def __setitem__(self, key, value):
        """Set data, automatically handling circular buffer indexing for the first dimension."""
        if isinstance(key, int):
            # Single event index - handle circular buffer logic
            # self._buffer._check_and_flush_if_needed(key)
            buffer_position = key % self._buffer.buffer_size
            self._array[buffer_position] = value
            self._buffer._number_written[self._buffer._key] += 1
            self._buffer._total_events_processed = max(self._buffer._total_events_processed, key + 1)
            # Track that this event position has been written to
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
                self._buffer._number_written[self._buffer._key] += 1
                self._buffer._total_events_processed = max(self._buffer._total_events_processed, first_idx + 1)
                # Track that this event position has been written to
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

    def __init__(self, buffer_size: int = 1000, flush_handler: Optional[BufferFlushHandler] = None):
        """
        Initialize the buffer.

        Args:
            buffer_size: Maximum number of events to store before flushing
            flush_handler: Handler for flush operations (optional)
        """
        self.buffer_size = buffer_size
        self.nevents = -1
        self.flush_handler = flush_handler

        # Internal storage
        self._buffer_arrays: Dict[str, BufferArray] = {}
        self._array_specs: Dict[str, Tuple] = {}  # Store (shape, dtype) for each array

        # Track buffer state - much simpler approach
        self._current_size = 0  # Number of events currently in buffer
        self._total_events_processed = 0
        self._buffer_start_event = 0
        self._current_event_positions = set()  # Track which event positions have been written to

        self._number_written: Dict[str, int] = {}
        self._key = None

    def _initialize_array(self, key: str, shape: Tuple, dtype: np.dtype):
        """Initialize a new array in the buffer with the given specifications."""
        full_shape = (self.buffer_size,) + shape[1:]  # Replace first dim with buffer_size
        array = np.zeros(full_shape, dtype=dtype)
        self._buffer_arrays[key] = BufferArray(array, self)
        self._array_specs[key] = (shape, dtype)
        self._number_written[key] = 0

        # if(self.nevents > -1):
        #     self.nevents = shape[0]
        # else:
        #     print('self.nevents = {}, shape[0] = {}'.format(self.nevents,shape[0]))
        #     assert self.nevents == shape[0]

    def _check_and_flush_if_needed(self, event_index: int):
        """Check if we need to flush before processing this event index."""
        # buffer_position = event_index % self.buffer_size

        do_flush = True
        for key,val in self._number_written.items():
            # print('\t-> {}, {}'.format(key,val))
            if(val != self.buffer_size):
                do_flush = False
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
        self._number_written[key] += copy_size

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
        if key not in self._buffer_arrays:
            arr_shape = (1,) + shape
            self._initialize_array(key, arr_shape, dtype)
        return self._buffer_arrays[key]

# Example usage demonstrating improved flushing behavior
if __name__ == "__main__":
    # Create a buffer with file flush handler
    flush_handler = DummyFlushHandler("physics_data.h5")
    buffer = Buffer(buffer_size=3, flush_handler=flush_handler)

    # Pre-create arrays with known maximum dimensions
    buffer.create_array('jets.N', dtype=np.int32)  # Scalar per event
    buffer.create_array('jets.Pmu', (10, 4), np.float64)  # Up to 10 jets, 4-momentum each
    buffer.create_array('event.weight', dtype=np.float64)  # Event weight

    print("Processing events with improved flushing logic...")
    print("Buffer will only flush when about to overwrite unflushed data\n")

    # Process events
    for event_index in range(8):
        # Simulate jet data for this event
        n_jets = np.random.randint(2, 6)
        jet_vectors = np.random.randn(n_jets, 4)  # (N_jets, 4-momentum)
        jet_ordering = np.arange(n_jets)  # Dummy ordering
        event_weight = np.random.uniform(0.5, 2.0)

        print(f"Processing event {event_index} (buffer pos {event_index % buffer.buffer_size})")

        # Fill all branches for this event
        buffer['jets.N'][event_index] = len(jet_vectors)

        embed_array_inplace(np.vstack([jet_vectors[i] for i in jet_ordering]),
                           buffer['jets.Pmu'][event_index])

        buffer['event.weight'][event_index] = event_weight

        print(f"  -> Jets: {n_jets}, buffer current size: {buffer._current_size}")
        print()

    # Final flush for remaining data
    print(f"Final flush of remaining {buffer._current_size} events:")
    buffer.flush()
    print("\nBuffer info:", buffer.get_buffer_info())