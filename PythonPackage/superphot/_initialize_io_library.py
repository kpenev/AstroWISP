"""Define a function for intefracing with the superphotio library."""

from ctypes import cdll, c_void_p, c_char_p
from ctypes.util import find_library

#Naming convention imitates the one by ctypes.
#pylint: disable=invalid-name

#Type checking place holders require no content.
#pylint: disable=too-few-public-methods
class _c_h5io_data_tree_p(c_void_p):
    """Placeholder for pointer to H5IODataTree opaque struct."""

#pylint: enable=invalid-name
#pylint: enable=too-few-public-methods

def _initialize_io_library():
    """Load and initalize the superphotio librray."""

    library_fname = find_library('superphotio')
    if library_fname is None:
        raise OSError('Unable to find the SuperPhot I/O library.')
    library = cdll.LoadLibrary(library_fname)

    library.create_result_tree.argtypes = [c_void_p, c_char_p]
    library.create_result_tree.restype = _c_h5io_data_tree_p

    library.destroy_result_tree.argtypes = [
        library.create_result_tree.restype
    ]
    library.destroy_result_tree.restype = None

    return library
