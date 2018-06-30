"""Define a function for intefracing with the superphotio library."""

from ctypes import\
    cdll,\
    c_void_p,\
    c_char_p,\
    c_bool,\
    c_double,\
    c_int,\
    c_short,\
    c_long,\
    c_byte,\
    c_uint,\
    c_ulong,\
    c_ushort,\
    c_ubyte,\
    pointer,\
    cast
from ctypes.util import find_library
import numpy

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

    library.query_result_tree.argtypes = [
        library.create_result_tree.restype,
        c_char_p,
        c_char_p,
        c_void_p
    ]
    library.query_result_tree.restype = c_bool

    library.free.argtypes = [c_void_p]
    library.free.restype = None

    return library

#Sufficient functionality to justify a class.
#pylint: disable=too-few-public-methods
class SuperPhotIOTree:
    """Interface for extracting entries from an IO::H5IODataTree."""

    library = _initialize_io_library()

    type_string = {c_bool: 'bool',
                   c_double: 'double',
                   c_int: 'int',
                   c_short: 'short',
                   c_long: 'long',
                   c_byte: 'char',
                   c_uint: 'uint',
                   c_ulong: 'ulong',
                   c_ushort: 'ushort',
                   c_ubyte: 'uchar'}

    def get(self, quantity, dtype=c_double, shape=(1,)):
        """
        Return the given quantity as a proper python object.

        Args:
            quantity (str):    The quantity to extract from the tree.

            dtype:    The data type of individual values of the quantity.

            shape (tulpe of ints):    The shape of the array of values
                to expect.

        Returns:
            numpy.ndarray(shape=shape, dtype=dtype):
                The values of the quantity. The return type is always an array,
                even for sintgle valued quantities. In the latter case, the
                shape is (1,).
        """

        print('Reading result quantity: ' + repr(quantity)
              +
              ', type: ' + repr(dtype)
              +
              ', shape: ' + repr(shape))
        byte_quantity = (quantity if isinstance(quantity, bytes)
                         else quantity.encode('ascii'))

        if dtype == str:
            library_result = pointer(c_char_p())
            defined = self.library.query_result_tree(
                self.library_tree,
                byte_quantity,
                b'str',
                cast(library_result, c_void_p)
            )
            result = library_result.contents.value.decode()
            print('Freeing library result')
            self.library.free(library_result.contents)
        else:
            result = numpy.empty(shape=shape, dtype=dtype)

            type_string_arg = self.type_string[dtype]
            for dim_size in shape:
                if dim_size != 1:
                    type_string_arg = '[' + type_string_arg + ']'
                    break

            defined = self.library.query_result_tree(
                self.library_tree,
                byte_quantity,
                type_string_arg.encode('ascii'),
                result.ctypes.data_as(c_void_p)
            )
        if not defined:
            raise KeyError(
                'Given result tree does not contain a quantity named: '
                +
                repr(quantity)
            )
        return result

    def __init__(self, library_configuration, version_info=''):
        """
        Create a tree with just the given configuration.

        Args:
            library_configuration:    The configuration object created by the
                library for the tool which will be using the configuration tree.

            version_info:    Information about the version of the
                tool/scripts/... using this tree. It is safe to leave this
                empty, if it is not required as an entry in the tree.

        Returns:
            None
        """

        self.library_tree = self.library.create_result_tree(
            library_configuration,
            (
                version_info if isinstance(version_info, bytes)
                else version_info.encode('ascii')
            )
        )

    def __del__(self):
        """Destroy the tree allocated by __init__."""

        print('Destroying result tree.')
        self.library.destroy_result_tree(self.library_tree)
#pylint: enable=too-few-public-methods
