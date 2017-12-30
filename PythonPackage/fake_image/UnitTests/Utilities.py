import unittest
from sys import float_info

class FloatTestCase(unittest.TestCase) :
    """Test case involving floating point comparisons."""

    tolerance = 10.0

    def approx_equal(self, v1, v2) :
        """Return True iff v1 ~= v2 to within tolerance."""

        return (abs(v1 - v2)
                <=
                self.tolerance * abs(v1 + v2) * float_info.epsilon)

    def set_tolerance(self, tolerance) :
        """
        Set the tolerance for floating point comparisons.

        Args:
            - tolerance:
                The maximum fractional difference in units of floating point
                epsilon to allow when comparing two floating point numbers.
                See assertApprox for more details.

        Returns: None
        """

        self.tolerance = tolerance

    def assertApprox(self, v1, v2, message = '') :
        """
        Assert that the two values are equal to wihin the current tolerance.

        The exact definition is: 

        |v1 - v2| <= tolerange * |v1 + v2|  * float_info.epsilon

        Args:
            - v1:
                The first of the two values to compare.
            - v2:
                The second of the two values to compare.

        Returns:
            True iff the two values are within tolerance of each other.
        """

        self.assertTrue(
            self._approx_equal(v1, v2),
            '%f !~ %f: %s' % (v1, v2, message)
        )
