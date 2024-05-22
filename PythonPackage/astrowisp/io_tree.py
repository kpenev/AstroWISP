"""Define a function for intefracing with the astrowisp shared library."""

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

from astrowisp._initialize_library import get_astrowisp_library

#Sufficient functionality to justify a class.
#pylint: disable=too-few-public-methods
class IOTree:
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

    def __init__(self, tool_or_configuration, version_info=''):
        """
        Create a tree with just the given configuration.

        Args:
            tool_or_configuration:    The configuration object created by the
                library for the tool which will be using the configuration tree
                or the tool itself.

            version_info:    Information about the version of the
                tool/scripts/... using this tree. It is safe to leave this
                empty, if it is not required as an entry in the tree.

        Returns:
            None
        """

        self._astrowisp_library = get_astrowisp_library()
        self.library_tree = self._astrowisp_library.create_result_tree(
            getattr(
                tool_or_configuration,
                '_library_configuration',
                tool_or_configuration
            ),
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
        num_quantities = self._astrowisp_library.list_tree_quantities(
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
            defined = self._astrowisp_library.query_result_tree(
                self.library_tree,
                byte_quantity,
                b'str',
                cast(library_result, c_void_p)
            )
            result = library_result.contents.value
            if result is None:
                defined = False
            else:
                result = result.decode()
            self._astrowisp_library.export_free(library_result.contents)
        else:
            result = numpy.empty(shape=shape, dtype=dtype)

            if dtype in self.type_string:
                type_string_arg = self.type_string[dtype]
            else:
                type_string_arg = dtype.str.lstrip('|')
            if shape is not None:
                type_string_arg = '[' + type_string_arg + ']'

            defined = self._astrowisp_library.query_result_tree(
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
        self._astrowisp_library.get_psf_map_variables(self.library_tree,
                                                      image_index,
                                                      result)
        return result

    def set_star_shape_map(self, grid, coefficients):
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

        self._astrowisp_library.update_result_tree(
            b'psffit.grid',
            (c_char_p * 1)(c_char_p(get_grid_str())),
            b'str',
            1,
            self.library_tree
        )
        c_coefficients = coefficients.astype(c_double, 'C')
        if coefficients.size:
            self._astrowisp_library.update_result_tree(
                b'psffit.psfmap',
                c_coefficients.ctypes.data_as(c_void_p),
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
                                       magnitude_1adu=None,
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

            star_shape_grid(2-D iterable):    The grid boundaries on which the
                star shape is being modeled.

            star_shape_map_terms(2-D numpy array):    The values of the terms
                required to evaluate the PSF/PRF map for each source. First
                dimension should iterate over sources and second over
                expansion terms.

            magnitude_1adu(float):    The magnitude that corresponds to a flux
                of 1ADU. Only required if relying on magnitudes to get star
                shape amplitudes.

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
            (
                'mag' in source_var_set
                and
                'mag_err' in source_var_set
                and
                magnitude_1adu is not None
            )
        )
        for prefix, prefix_vars in [('projsrc', ['x',
                                                 'y',
                                                 'enabled']),
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
                    assert var_name == 'enabled' or dtype.kind == 'f'
                    assert var_name == 'enabled' or dtype.itemsize == 8
                    c_data = source_data[var_name].astype(
                        c_double,
                        order='C'
                    )
                    self._astrowisp_library.update_result_tree(
                        '.'.join(
                            [prefix, var_name, image_index_str]
                        ).encode('ascii'),
                        c_data.ctypes.data_as(
                            c_void_p
                        ),
                        b'double',
                        source_data.shape[0],
                        self.library_tree
                    )

        image_index_str = image_index_str.encode('ascii')
        self._astrowisp_library.update_result_tree(
            b'projsrc.srcid.name.' + image_index_str,
            (c_char_p * source_data.shape[0])(*source_data['ID']),
            b'str',
            source_data.shape[0],
            self.library_tree
        )

        c_data = source_data['bg'].astype(c_double)
        self._astrowisp_library.update_result_tree(
            b'bg.value.' + image_index_str,
            c_data.ctypes.data_as(c_void_p),
            b'double',
            source_data.shape[0],
            self.library_tree
        )
        c_data = source_data['bg_err'].astype(c_double)
        self._astrowisp_library.update_result_tree(
            b'bg.error.' + image_index_str,
            c_data.ctypes.data_as(c_void_p),
            b'double',
            source_data.shape[0],
            self.library_tree
        )
        c_data = source_data['bg_npix'].astype(c_uint)
        self._astrowisp_library.update_result_tree(
            b'bg.npix.' + image_index_str,
            c_data.ctypes.data_as(c_void_p),
            b'uint',
            source_data.shape[0],
            self.library_tree
        )
        self._astrowisp_library.update_result_tree(
            b'psffit.magnitude_1adu',
            numpy.array([magnitude_1adu],
                        dtype=c_double).ctypes.data_as(c_void_p),
            b'double',
            1,
            self.library_tree
        )
        self._astrowisp_library.update_result_tree(
            b'psffit.model',
            (c_char_p * 1)(b'bicubic'),
            b'str',
            1,
            self.library_tree
        )
        assert star_shape_map_terms.shape == (
            source_data.size,
            star_shape_map_coefficients.shape[-1]
        )
        c_data = star_shape_map_terms.astype(c_double)
        self._astrowisp_library.update_result_tree(
            b'psffit.terms.' + image_index_str,
            c_data.ctypes.data_as(c_void_p),
            b'double',
            star_shape_map_terms.size,
            self.library_tree
        )

        self.set_star_shape_map(star_shape_grid, star_shape_map_coefficients)

    def __del__(self):
        """Destroy the tree allocated by __init__."""

        self._astrowisp_library.destroy_result_tree(self.library_tree)
#pylint: enable=too-few-public-methods
