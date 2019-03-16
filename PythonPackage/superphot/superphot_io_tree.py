"""Define a function for intefracing with the superphotio library."""

from ctypes import\
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
    c_char_p,\
    c_void_p,\
    pointer,\
    POINTER,\
    byref,\
    cast

import numpy

from superphot._initialize_library import superphot_library

#Sufficient functionality to justify a class.
#pylint: disable=too-few-public-methods
class SuperPhotIOTree:
    """Interface for extracting entries from an IO::H5IODataTree."""

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

        self.library_tree = superphot_library.create_result_tree(
            library_configuration,
            (
                version_info if isinstance(version_info, bytes)
                else version_info.encode('ascii')
            )
        )

    def defined_quantity_names(self):
        """
        Return a list of the quantities with non-empty values in the tree.

        Args:
            None

        Returns:
            [str]:
                The full names of the quantities with available data using dot
                as a separator between leves in the tree.
        """

        library_type = POINTER(c_char_p)
        library_result = library_type()
        num_quantities = superphot_library.list_tree_quantities(
            self.library_tree,
            byref(library_result)
        )
        return [library_result[index].decode()
                for index in range(num_quantities)]

    def get(self, quantity, dtype=c_double, shape=None):
        """
        Return the given quantity as a proper python object.

        Args:
            quantity (str):    The quantity to extract from the tree.

            dtype:    The data type of individual values of the quantity.

            shape (tuple of ints):    The shape of the array of values
                to expect. Use None for scalar quantities.

        Returns:
            numpy.ndarray(shape=shape, dtype=dtype):
                The values of the quantity. The return type is always an array,
                even for sintgle valued quantities. In the latter case, the
                shape is (1,).
        """

        byte_quantity = (quantity if isinstance(quantity, bytes)
                         else quantity.encode('ascii'))

        if dtype == str:
            library_result = pointer(c_char_p())
            defined = superphot_library.query_result_tree(
                self.library_tree,
                byte_quantity,
                b'str',
                cast(library_result, c_void_p)
            )
            result = library_result.contents.value.decode()
            superphot_library.free(library_result.contents)
        else:
            result = numpy.empty(shape=shape, dtype=dtype)

            type_string_arg = self.type_string[dtype]
            if shape is not None:
                type_string_arg = '[' + type_string_arg + ']'

            defined = superphot_library.query_result_tree(
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

    def get_psfmap_variables(self, image_index, num_variables, num_sources):
        """
        Return the values of the PSF map variables for all sources in an image.

        Args:
            image_index:    The index of the image for which to return the
                values of the variables as supplied to PSF fitting.

            num_variables:    The number of variables used for PSF fitting.

            num_sources:    The number of sources in the selected image.

        Returns:
            numpy.ndarray(dtype=float, shape=(num_variables, num_sources):
                Array with records named as the PSF map variables and entries
                containing the values of the variables for all sources in the
                image identified by image_index.
        """

        result = numpy.empty(dtype=float, shape=(num_variables, num_sources))
        superphot_library.get_psf_map_variables(self.library_tree,
                                                image_index,
                                                result)
        return result

    def __del__(self):
        """Destroy the tree allocated by __init__."""

        superphot_library.destroy_result_tree(self.library_tree)
#pylint: enable=too-few-public-methods
