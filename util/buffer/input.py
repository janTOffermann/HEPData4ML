
import h5py as h5
import uproot as ur
import awkward as ak
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
        self.nevents = dataset.shape[0]
        self.dtype = dataset.dtype
        self.ndim = dataset.ndim

        # Buffer state
        self._buffer = None
        self._buffer_start = -1
        self._buffer_end = -1

    def _load_buffer(self, event_idx: int) -> None:
        # Calculate buffer range
        buffer_start = (event_idx // self.buffer_size) * self.buffer_size
        buffer_end = min(buffer_start + self.buffer_size, self.nevents)

        # Only reload if we need a different buffer range
        if buffer_start != self._buffer_start:
            self._buffer = self.dataset[buffer_start:buffer_end]
            self._buffer_start = buffer_start
            self._buffer_end = buffer_end

    def __getitem__(self, key):
        if isinstance(key, int):
            # Single event access
            if key < 0:
                key = self.nevents + key  # Handle negative indexing

            if key < 0 or key >= self.nevents:
                raise IndexError(f"Index {key} out of bounds for axis 0 with size {self.nevents}")

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
                    event_idx = self.nevents + event_idx

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
        return self.nevents

    def __repr__(self):
        return f"BufferedHDF5Array(nevents={self.nevents}, dtype={self.dtype}, buffer_size={self.buffer_size})"

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

class UprootTreeLoader:
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


####################
# Experimental
####################

class LazyBranch:
    """Lazy accessor for a single branch with batch caching."""

    def __init__(self, loader, branch_name):
        self._loader = loader
        self._branch_name = branch_name

    def __getitem__(self, key):
        """Access entries from this branch, loading batches as needed."""
        if isinstance(key, int):
            batch, local_idx = self._loader._get_batch_for_index(key)
            return batch[self._branch_name][local_idx]

        elif isinstance(key, slice):
            # For slices, delegate to loader's slice logic
            sliced_data = self._loader[key]
            return sliced_data[self._branch_name]

        else:
            raise TypeError(f"Invalid key type: {type(key)}")

    def __len__(self):
        return len(self._loader)


class UprootBatchLoader:
    def __init__(self, files, tree_path, expressions, batch_size=100000):
        """
        Initialize a batched ROOT file loader using uproot.

        Parameters
        ----------
        files : str or list of str
            Path(s) to ROOT file(s)
        tree_path : str
            Path to TTree within the ROOT file(s)
        expressions : list of str
            Branch names to load
        batch_size : int, optional
            Number of entries to load per batch (default: 100000)
        """
        self.files = files
        self.tree_path = tree_path
        self.batch_size = batch_size

        # Get field info by checking what actually exists
        self.fields = None
        with ur.open(files[0]) as f:
            self.fields = [expr for expr in expressions if expr in f[tree_path]]

        # Initialize iteration state
        self._iterator = None
        self._current_batch = None
        self._current_batch_start = 0
        self._total_entries = None

        # Create the iterator
        self._reset_iterator()

    def _reset_iterator(self):
        """Reset the batch iterator to the beginning."""
        self._iterator = ur.iterate(
            self.files,
            expressions=self.fields,
            tree_path=self.tree_path,
            library="ak",
            step_size=self.batch_size
        )
        self._current_batch = None
        self._current_batch_start = 0

    def _get_batch_for_index(self, global_idx):
        """
        Load the batch containing the given global index.

        Parameters
        ----------
        global_idx : int
            Global index across all files

        Returns
        -------
        tuple
            (batch_data, local_index_in_batch)
        """
        # If we need to reset (going backwards or past current batch)
        if self._current_batch is None or global_idx < self._current_batch_start:
            self._reset_iterator()

        # Iterate until we find the right batch
        while True:
            batch_end = self._current_batch_start + (
                len(self._current_batch) if self._current_batch is not None else 0
            )

            if self._current_batch is not None and global_idx < batch_end:
                # Found the right batch
                local_idx = global_idx - self._current_batch_start
                return self._current_batch, local_idx

            # Load next batch
            try:
                self._current_batch_start = batch_end
                self._current_batch = next(self._iterator)
            except StopIteration:
                raise IndexError(f"Index {global_idx} out of range")

    def __getitem__(self, key):
        """
        Access branches or entries.

        Parameters
        ----------
        key : str or int or slice
            If str: return LazyBranch accessor for that branch
            If int: return entry at that index across all branches
            If slice: return slice of entries across all branches

        Returns
        -------
        LazyBranch or awkward.Array
            Requested data
        """
        if isinstance(key, str):
            # Return a lazy branch accessor
            return LazyBranch(self, key)

        elif isinstance(key, int):
            # Single entry access - return all branches for this event
            batch, local_idx = self._get_batch_for_index(key)
            return batch[local_idx]

        elif isinstance(key, slice):
            # Slice access - load relevant batches
            start, stop, step = key.indices(len(self))

            # Reset and collect relevant batches
            self._reset_iterator()
            result_batches = []
            current_pos = 0

            for batch in self._iterator:
                batch_end = current_pos + len(batch)

                # Check if this batch overlaps with our slice
                if batch_end > start and current_pos < stop:
                    # Calculate local slice within this batch
                    local_start = max(0, start - current_pos)
                    local_stop = min(len(batch), stop - current_pos)
                    result_batches.append(batch[local_start:local_stop])

                current_pos = batch_end

                if current_pos >= stop:
                    break

            result = ak.concatenate(result_batches) if result_batches else ak.Array([])
            return result[::step] if step != 1 else result

        else:
            raise TypeError(f"Invalid key type: {type(key)}")

    def __len__(self):
        """Return total number of entries across all files."""
        if self._total_entries is None:
            self._total_entries = 0
            for file_path in (self.files if isinstance(self.files, list) else [self.files]):
                with ur.open(file_path) as f:
                    self._total_entries += f[self.tree_path].num_entries
        return self._total_entries

    def iterate_batches(self):
        """
        Iterate over batches explicitly.

        Yields
        ------
        awkward.Array
            Each batch of data

        Examples
        --------
        >>> loader = UprootBatchLoader(files, "Delphes", ["Jet.PT", "Jet.Eta"])
        >>> for batch in loader.iterate_batches():
        ...     print(f"Processing {len(batch)} events")
        ...     process_jets(batch["Jet.PT"], batch["Jet.Eta"])
        """
        self._reset_iterator()
        for batch in self._iterator:
            yield batch