"""Define a class for imposing smooth depndence of one quantity on others."""

from ctypes import\
    cdll,\
    c_double,\
    c_char_p,\
    c_uint,\
    POINTER,\
    byref
from ctypes.util import find_library
import numpy.ctypeslib

def _initialize_library():
    """Prepare the SuperPhot psf library for use."""

    psf_library_fname = find_library('superphotpsf')
    io_library_fname = find_library('superphotio')
    if psf_library_fname is None:
        raise OSError('Unable to find the SuperPhot PSF library.')
    if io_library_fname is None:
        raise OSError('Unable to find the SuperPhot I/O library.')
    cdll.LoadLibrary(io_library_fname)
    library = cdll.LoadLibrary(psf_library_fname)

    library.expand_term_expression.argtypes = [
        c_char_p,
        POINTER(POINTER(c_char_p)),
        POINTER(c_uint)
    ]
    library.expand_term_expression.restype = None

    library.free_term_list.argtypes = [POINTER(c_char_p), c_uint]
    library.free_term_list.restype = None

    library.evaluate_terms.argtypes = [
        POINTER(c_char_p),
        c_uint,
        POINTER(c_char_p),
        POINTER(POINTER(c_double)),
        c_uint,
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=2,
                                  flags='C_CONTIGUOUS')
    ]
    return library

class SmoothDependence:
    r"""
    Utilities for imposing smooth depnedence of one quantity on others.

    Attributes:
        expression (str):    The expression defining the dependence of
            parameters on a set of variables.

        terms ([str]):    The terms the depnedence expression expands to.

        parameters (dict-like, record array or #-D array, #>=2):    The
            parameters that depend smoothly on the variables. If dict-like or
            record array, indexed by the names of the parameters and containing
            the corresponding expansion coefficients. If array, the first
            dimension is assumed to contain the expansion coefficients.

    Examples:

        * Find all terms combining up to x^3 and y^4, with only even powers of y
          allowed:

            >>> from superphot import SmoothDependence
            >>> term_list = SmoothDependence.expand_expression('O3{x} * O2{y^2}')

        * Evaluate all terms in a polynomial of x and y of up to combined third
          order:

          >>> from superphot import SmoothDependence
          >>> import numpy
          >>> term_values = SmoothDependence.evaluate_terms(
          >>>     'O3{x, y}',
          >>>     x=numpy.arange(10.0),
          >>>     y=numpy.arange(-10.0, 0.0)
          >>> )

        * Define two unnamed parameters smoothly depending on x:

              * :math:`2.0 x + 4.0 x^2 + 6.0 x^3`
              * :math:`1.0 + 3.0 x + 5.0 x^2 + 7.0 x^3`

          and evaluate them for :math:`x \in \{0, 0.1, 0.2, 0.3, 0.4\}`::

          >>> from superphot import SmoothDependence
          >>> import numpy
          >>> parameters = numpy.arange(8.0).reshape((4, 2))
          >>> smooth = SmoothDependence('O3{x}', parameters)
          >>> param_values = smooth(x=numpy.arange(0.0, 0.45, 0.1))

        * Same as above, but this time name the parameters `a` and `b`::

          >>> from superphot import SmoothDependence
          >>> import numpy
          >>> parameters = numpy.empty(4, dtype = [('a', float), ('b', float)])
          >>> parameters['a'] = numpy.arange(0.0, 7.5, 2.0)
          >>> parameters['b'] = numpy.arange(1.0, 7.5, 2.0)
          >>> smooth = SmoothDependence('O3{x}', parameters)
          >>> param_values = smooth(x=numpy.arange(0.0, 0.45, 0.1))
    """

    _library = _initialize_library()

    @classmethod
    def expand_expression(cls, expression):
        """
        Return the terms expression expands to per SuperPhot smooth grammar.

        Args:
            expression(str or bytes):    The expression to parse.

        Returns:
            [str]:
                List of the term the expression expands to.
        """

        terms_type = POINTER(c_char_p)
        terms = terms_type()
        num_terms = c_uint()
        c_expression = c_char_p(expression if isinstance(expression, bytes)
                                else expression.encode('ascii'))
        cls._library.expand_term_expression(c_expression,
                                            byref(terms),
                                            byref(num_terms))
        result = [terms[term_index].decode()
                  for term_index in range(num_terms.value)]

        cls._library.free_term_list(terms, num_terms)
        return result

    @classmethod
    def evaluate_terms(cls, terms, **indep_var):
        """
        Evaluate the given terms or a term expression for the given sources.

        Args:
            terms (str or [str]):    Either a single string giving an expression
                which is expanded to a list of terms using expand_expression(),
                or directly a list of terms to evaluate.

            **indep_var (1D numpy array):    The values of the independent
                variables involved in terms. The argument names must match the
                names used in `terms`.

        Returns:
            2D numpy array:
                The values of the terms for each entry in indep_var. Each column
                contains a different term.
        """

        if isinstance(terms, (str, bytes)):
            terms = cls.expand_expression(terms)
        term_list = (c_char_p * len(terms))(*(t.encode('ascii') for t in terms))

        num_variables = len(indep_var)
        num_terms = len(terms)
        num_sources = len(next(iter(indep_var.values())))

        #Only need names not values.
        #pylint: disable=consider-iterating-dictionary
        variable_names = (c_char_p * num_variables)(
            *(v.encode('ascii') for v in indep_var.keys())
        )
        #pylint: enable=consider-iterating-dictionary

        contiguous_vars = [numpy.require(values, c_double, ('C_CONTIGUOUS',))
                           for values in indep_var.values()]
        variable_values = (POINTER(c_double) * num_variables)(
            *(var.ctypes.data_as(POINTER(c_double)) for var in contiguous_vars)
        )

        result = numpy.empty((num_sources, num_terms), dtype=c_double)

        cls._library.evaluate_terms(term_list,
                                    num_terms,
                                    variable_names,
                                    variable_values,
                                    num_variables,
                                    num_sources,
                                    result)
        return result

    def __init__(self, expression, parameters):
        """
        Set-up a smooth depnednece of a set of parameters on a set of variables.

        Args:
            expression:    See :attr:`expression` attribute.

            parameters:    See :attr:`parameters` attribute.

        Returns:
            None
        """

        self.expression = expression
        self.terms = self.expand_expression(expression)
        self.parameters = parameters

        assert parameters.shape[0] == len(self.terms)
        assert (
            (parameters.ndim >= 2 and parameters.dtype.names is None)
            or
            (parameters.ndim == 1 and parameters.dtype.names is not None)
        )

    def __call__(self, **indep_var):
        """
        Calculate the :attr:`parameters` for the given variable values.

        Args:
            **indep_var (1D arrays):    The values of the independent variables
                to evaluate the dependence at. Must all have the same length.

        Returns:
            (same type as :attr:`parameters` attribute):
                The values of the parameters evaluated for the indepnedent
                variable values supplied.
        """

        num_sources = len(next(iter(indep_var.values())))
        num_terms = len(self.terms)

        result = numpy.empty(shape=(num_sources,) + self.parameters.shape[1:],
                             dtype=self.parameters.dtype)

        if result.dtype.names:
            num_parameters = len(self.parameters.dtype.names)

            #pylint false positive
            #pylint: disable=unused-variable
            write_to = result.view(float).reshape(
                (num_sources, num_parameters)
            )
            #pylint: enable=unused-variable

            parameter_view = self.parameters.view(float).reshape(
                (num_terms, num_parameters)
            )
        else:
            write_to = result
            parameter_view = self.parameters

        term_values = self.evaluate_terms(self.terms, **indep_var)
        write_to[:] = numpy.tensordot(
            term_values,
            parameter_view,
            1
        )

        return result

if __name__ == '__main__':
    test_terms = SmoothDependence.expand_expression('O3{x} * O2{y^2}')
    print(test_terms)
    print(
        repr(
            SmoothDependence.evaluate_terms('O3{x, y}',
                                            x=numpy.arange(10.0),
                                            y=numpy.arange(-10.0, 0.0))
        )
    )

    named_param = numpy.empty(4, dtype=[('a', float), ('b', float)])
    named_param['a'] = numpy.arange(0.0, 7.5, 2.0)
    named_param['b'] = numpy.arange(1.0, 7.5, 2.0)
    for param in [numpy.arange(8.0).reshape((4, 2)),
                  named_param]:
        xcubed = SmoothDependence('O3{x}', param)
        print(repr(xcubed(x=numpy.arange(0.0, 0.45, 0.1))))
