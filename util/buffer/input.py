
import h5py as h5
import uproot as ur
from typing import List

# TODO: Consider making more general classes, and using inheritance.
#       Ultimately should support reading from ROOT files as well.

class BufferedHDF5Array:
    """
    A numpy-array-like interface for buffered reading of HDF5 datasets.
    """

    def __init__(self, dataset: h5.Dataset, buffer_size: int = 10000):

        self.dataset = dataset
        self.buffer_size = buffer_size
        self.shape = dataset.shape
        self.dtype = dataset.dtype
        self.ndim = dataset.ndim

        # Buffer state
        self._buffer = None
        self._buffer_start = -1
        self._buffer_end = -1

    def _load_buffer(self, event_idx: int) -> None:
        # Calculate buffer range
        buffer_start = (event_idx // self.buffer_size) * self.buffer_size
        buffer_end = min(buffer_start + self.buffer_size, self.shape[0])

        # Only reload if we need a different buffer range
        if buffer_start != self._buffer_start:
            self._buffer = self.dataset[buffer_start:buffer_end]
            self._buffer_start = buffer_start
            self._buffer_end = buffer_end

    def __getitem__(self, key):
        if isinstance(key, int):
            # Single event access
            if key < 0:
                key = self.shape[0] + key  # Handle negative indexing

            if key < 0 or key >= self.shape[0]:
                raise IndexError(f"Index {key} out of bounds for axis 0 with size {self.shape[0]}")

            self._load_buffer(key)
            buffer_idx = key - self._buffer_start
            return self._buffer[buffer_idx]

        elif isinstance(key, slice):
            # Slice access - for now, just read directly from dataset
            return self.dataset[key]

        elif isinstance(key, tuple):
            # Multi-dimensional indexing
            event_idx = key[0]
            remaining_key = key[1:]

            if isinstance(event_idx, int):
                if event_idx < 0:
                    event_idx = self.shape[0] + event_idx

                self._load_buffer(event_idx)
                buffer_idx = event_idx - self._buffer_start
                return self._buffer[buffer_idx][remaining_key]
            else:
                # Complex indexing, fall back to direct dataset access
                return self.dataset[key]
        else:
            # Other indexing types, fall back to dataset
            return self.dataset[key]

    def __len__(self):
        return self.shape[0]

    def __repr__(self):
        return f"BufferedHDF5Array(shape={self.shape}, dtype={self.dtype}, buffer_size={self.buffer_size})"

class BufferedHDF5Reader:
    """
    A class for efficiently reading multiple HDF5 datasets with coordinated buffering.
    """

    def __init__(self, h5_file: str, collections: List[str], buffer_size: int = 10000):

        self.h5_file = h5_file
        self.collections = collections
        self.buffer_size = buffer_size

        # Open file and create buffered arrays
        self._file = h5.File(h5_file, 'r')
        self._arrays = {}

        # Get number of events from first collection
        self.nevents = self._file[collections[0]].shape[0]

        # Verify all collections have same number of events
        for key in collections:
            if self._file[key].shape[0] != self.nevents:
                raise ValueError(f"Collection '{key}' has {self._file[key].shape[0]} events, "
                               f"expected {self.nevents}")

        # Create buffered arrays
        for key in collections:
            self._arrays[key] = BufferedHDF5Array(self._file[key], buffer_size)

    def __getitem__(self, key: str) -> BufferedHDF5Array:
        return self._arrays[key]

    def __contains__(self, key: str) -> bool:
        return key in self._arrays

    def keys(self):
        return self._arrays.keys()

    def items(self):
        return self._arrays.items()

    def values(self):
        return self._arrays.values()

    def close(self):
        if hasattr(self, '_file'):
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _ensure_event_buffered(self, event_idx: int):
        """
        Internal method to ensure all arrays have the specified event buffered.

        This is called automatically by BufferedHDF5Array.__getitem__() to
        coordinate buffering across all collections.
        """
        # Check if any array needs a buffer reload for this event
        arrays_needing_reload = []
        for key, array in self._arrays.items():
            if array._needs_buffer_reload(event_idx):
                arrays_needing_reload.append((key, array))

        # If any need reloading, reload all of them to stay coordinated
        if arrays_needing_reload:
            for key, array in arrays_needing_reload:
                array._load_buffer(event_idx)

########################################################
# Uproot/dask-related stuff, for use with Delphes reading.
########################################################

class UprootBatchLoader:
    def __init__(self, files, tree_path, expressions):
        self.files = files
        self.tree_path = tree_path

        # Get field info by checking what actually exists
        self.fields = None
        with ur.open(files[0]) as f:
            self.fields = [expr for expr in expressions if expr in f[tree_path]]

        # Load all branches at once
        self._data = ur.concatenate(
            files,
            expressions=self.fields,
            tree_path=tree_path,
            library="ak"
        )

    def __getitem__(self, branch_name):
        return self._data[branch_name]

    def __len__(self):
        if self._data is not None:
            return len(self._data)
        else:
            total = 0
            for file_path in self.files:
                with ur.open(file_path) as f:
                    total += f[self.tree_path].num_entries
            return total
