"""A collection of general use testing utilities."""

import unittest
from sys import float_info


class FloatTestCase(unittest.TestCase):
    """Test case involving floating point comparisons."""

    relative_tolerance = 10.0
    absolute_tolerance = 0.0

    def approx_equal(self, value_1, value_2):
        """Return True iff value_1 ~= value_2 to within tolerance."""

        return abs(value_1 - value_2) <= (
            self.relative_tolerance
            * abs(value_1 + value_2)
            * float_info.epsilon
            + self.absolute_tolerance
        )

    def set_tolerance(self, relative, absolute=0.0):
        """
        Set the tolerance for floating point comparisons.

        Args:
            tolerance:    The maximum fractional difference in units of floating
                point epsilon to allow when comparing two floating point
                numbers. See assertApprox for more details.

        Returns:
            None
        """

        self.relative_tolerance = relative
        self.absolute_tolerance = absolute

    # Following standard unittest assert naming convections
    # pylint: disable=invalid-name
    def assertApprox(self, value_1, value_2, message=""):
        r"""
        Assert that the two values are equal to wihin the current tolerance.

        Notes:
            The exact definition is:

            \|value_1 - value_2\|
            <=
            tolerange * \|value_1 + value_2\|  * float_info.epsilon

        Args:
            value_1:    The first of the two values to compare.

            value_2:    The second of the two values to compare.

        Returns:
            True iff the two values are within tolerance of each other.
        """

        self.assertTrue(
            self.approx_equal(value_1, value_2),
            f"{value_1!r} !~ {value_2!r} (diff {value_1 - value_2}): {message}",
        )

    def assertApproxPandas(self, expected, testing, testing_what=""):
        """Assert that the two dataframes match columns, types and values."""

        if hasattr(self, '_logger'):
            self._logger.debug("Expected :\n%s", repr(expected))
            self._logger.debug("Got :\n%s", repr(testing))

        self.assertTrue(
            (
                expected.columns.size == testing.columns.size
                and (expected.columns == testing.columns).all()
            ),
            f"{testing_what}: column differ from expected:\n\t"
            f"{testing.columns!r}\n\tinstead of\n\t{expected.columns!r}",
        )
        for column in expected.columns:
            expected_col = expected[column]
            testing_col = testing[column]
            self.assertTrue(
                expected_col.dtype == testing_col.dtype,
                f"Column types mismatch: {testing_col.dtype!r} instead of "
                f"{expected_col.dtype!r}",
            )
            mismatch_message = (
                f"{testing_what}: Column {column!r} mismatch:"
                f"\n\t{testing_col!r}\n\tinstead of\n\t{expected_col!r}"
            )

            if expected_col.dtype.kind in ["i", "u"]:
                self.assertTrue(
                    expected_col.equals(testing_col), mismatch_message
                )
            else:
                assert expected_col.dtype.kind == "f"
            for expected_val, testing_val in zip(expected_col, testing_col):
                self.assertApprox(expected_val, testing_val, mismatch_message)

    # pylint: enable=invalid-name
