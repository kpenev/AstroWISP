"""Define a function for intefracing with the superphotio library."""

from ctypes import\
    c_bool,\
    c_double,\
    c_int,\
    c_short,\
    c_long,\
    c_byte,\
    c_ubyte,\
    c_char,\
    c_uint,\
    c_ulong,\
    c_ushort,\
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
                   c_ubyte: 'uchar',
                   c_char: 'char',
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

            if dtype in self.type_string:
                type_string_arg = self.type_string[dtype]
            else:
                type_string_arg = dtype.str.lstrip('|')
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

    def set_star_shape_map_variables(self, source_data, image_index):
        """
        Add the variables the star shape map depends on to the tree.

        Args:
            source_data (structured numpy.array):    The data to build the
                variables from. All floating point fields except 'bg' and
                'bg_err' are added as PSF map variables.

        Returns:
            None
        """

        def get_variable_names():
            """Identify and return the variable names directly as c_char_p."""

            result = []
            for var_name in source_data.dtype.names:
                if (
                        source_data[var_name].dtype.kind == 'f'
                        and
                        var_name not in ['bg', 'bg_err',
                                         'flux', 'flux_err',
                                         'mag', 'mag_err']
                ):
                    result.append(var_name)
            return result

        def get_variable_values(variable_names):
            """Return a properly laid out array of the variable values."""

            result = numpy.empty(shape=(len(variable_names), source_data.shape[0]),
                                 dtype=float)
            for var_index, var_name in enumerate(variable_names):
                result[var_index] = source_data[var_name]
            return result

        variable_names = get_variable_names()
        c_variable_names = (c_char_p * len(variable_names))(
            *(c_char_p(var_name.encode('ascii')) for var_name in variable_names)
        )
        superphot_library.set_psf_map_variables(
            c_variable_names,
            get_variable_values(variable_names),
            len(variable_names),
            source_data.shape[0],
            image_index,
            self.library_tree
        )

    def set_star_shape_map(self, grid, map_terms, coefficients):
        """
        Add to tree all entries that define the star shape map.

        Args:
            grid (2-D iterable):    The grid on which the star shape is
                represented.

            map_terms (str):    The expression defining the terms the star shape
                parameters depend on.

            coefficients (4-D numpy.array):    The coefficients of the map. See
                :class:fit_star_shape for details.

        Returns:
            None
        """

        def get_grid_str():
            """Return the ascii representation of the grid to add to self."""

            return ';'.join(
                [
                    ', '.join([repr(boundary) for boundary in sub_grid])
                    for sub_grid in grid
                ]
            ).encode('ascii')

        superphot_library.update_result_tree(
            b'psffit.grid',
            (c_char_p * 1)(c_char_p(get_grid_str())),
            b'str',
            1,
            self.library_tree
        )
        superphot_library.update_result_tree(
            b'psffit.terms',
            (c_char_p * 1)(
                c_char_p(
                    map_terms if isinstance(map_terms, bytes)
                    else map_terms.encode('ascii')
                )
            ),
            b'str',
            1,
            self.library_tree
        )
        superphot_library.update_result_tree(
            b'psffit.psfmap',
            coefficients.astype(c_double, 'C').ctypes.data_as(c_void_p),
            b'double',
            coefficients.size,
            self.library_tree
        )

    def set_aperture_photometry_inputs(self,
                                       *,
                                       source_data,
                                       star_shape_grid,
                                       star_shape_map_terms,
                                       star_shape_map_coefficients,
                                       image_index=0):
        """
        Add to the tree all the information required for aperture photometry.

        Args:
            source_data(structured numpy.array):    Should contain informaiton
                about all sources to do apreture photometry on as fields. At
                least the following floating point fields must be present: `x`,
                `y`, `bg`, `bg_err`, any variables used by the PSF map, and
                either `flux` and `flux_err` or `mag` and `mag_err`. It must
                also contain a string field `id` of source IDs and an unsigned
                integer field `bg_npix`.

            star_shape_grid:    The grid boundaries on which the star shape is
                being modeled.

            star_shape_map_terms:    The expression defining all terms to include in
                the star shape dependence.

            star_shape_map_coefficients(4-D numpy.array):    The coefficients
                in front of all terms. See bicubic PSf model for details.

        Returns:
            None
        """

        image_index_str = str(image_index)
        source_var_set = set(source_data.dtype.names)
        assert (
            ('flux' in source_var_set and 'flux_err' in source_var_set)
            or
            ('mag' in source_var_set and 'mag_err' in source_var_set)
        )
        for prefix, prefix_vars in [('projsrc', ['x',
                                                 'y']),
                                    ('psffit', ['flux',
                                                'flux_err',
                                                'mag',
                                                'mag_err'])]:
            for var_name in prefix_vars:
                if (
                        prefix == 'projsrc'
                        or
                        var_name in source_data.dtype.names
                ):
                    dtype = source_data[var_name].dtype
                    assert dtype.kind == 'f'
                    assert dtype.itemsize == 8
                    superphot_library.update_result_tree(
                        '.'.join(
                            [prefix, var_name, image_index_str]
                        ).encode('ascii'),
                        source_data[var_name].astype(
                            c_double
                        ).ctypes.data_as(
                            c_void_p
                        ),
                        b'double',
                        source_data.shape[0],
                        self.library_tree
                    )

        superphot_library.update_result_tree(
            b'projsrc.srcid.name.' + image_index_str.encode('ascii'),
            (c_char_p * source_data.shape[0])(*source_data['id']),
            b'str',
            source_data.shape[0],
            self.library_tree
        )

        superphot_library.update_result_tree(
            b'bg.value.' + image_index_str.encode('ascii'),
            source_data['bg'].astype(c_double).ctypes.data_as(c_void_p),
            b'double',
            source_data.shape[0],
            self.library_tree
        )
        superphot_library.update_result_tree(
            b'bg.error.' + image_index_str.encode('ascii'),
            source_data['bg_err'].astype(c_double).ctypes.data_as(c_void_p),
            b'double',
            source_data.shape[0],
            self.library_tree
        )
        superphot_library.update_result_tree(
            b'bg.npix.' + image_index_str.encode('ascii'),
            source_data['bg_npix'].astype(c_uint).ctypes.data_as(c_void_p),
            b'uint',
            source_data.shape[0],
            self.library_tree
        )

        self.set_star_shape_map_variables(source_data, image_index)

        self.set_star_shape_map(star_shape_grid,
                                star_shape_map_terms,
                                star_shape_map_coefficients)

    def __del__(self):
        """Destroy the tree allocated by __init__."""

        superphot_library.destroy_result_tree(self.library_tree)
#pylint: enable=too-few-public-methods
