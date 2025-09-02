# Note: Functions here use our custom VectorCalcs C++/ROOT library.
#       That library must be loaded, otherwise funcs will not work.
import numpy as np

def embed_array(array, target_shape,padding_value=0):
    """
    A generic function for embedding an input array into some target shape.
    The array will be truncated or padded as needed.

    """
    # Create the output array
    result = np.full(target_shape, padding_value, dtype=array.dtype)

    # Calculate the effective shape (minimum dimensions)
    effective_shape = tuple(min(s, t) for s, t in zip(array.shape, target_shape))

    slices = tuple(slice(0, dim) for dim in effective_shape)
    result[slices] = array[slices]
    return result

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