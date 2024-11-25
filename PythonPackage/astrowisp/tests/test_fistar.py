#!/usr/bin/env python3

"""Test source extraction using built-in fistar executable."""

from os import path
from subprocess import Popen, PIPE
from functools import partial

import unittest
from pandas import read_csv

from astrowisp import fistar_path
from astrowisp.utils.file_utilities import get_unpacked_fits
from astrowisp.tests.utilities import FloatTestCase

_test_data_dir = path.join(path.dirname(path.abspath(__file__)),
                           'test_data')

class TestFistar(FloatTestCase):
    """Test cases for the fistar executable."""

    def test_xo1(self):
        """Check if extracting sources from XO-1 image matches expected."""

        parse_result = partial(read_csv, sep=r'\s+', comment='#', header=None)

        expected = parse_result(
            path.join(_test_data_dir, 'XO1_test_img.fistar')
        )
        print(f'Expected:\n{expected!r}')
        with get_unpacked_fits(
            path.join(_test_data_dir, 'XO1_test_img.fits')
        ) as unpacked_fname:
            with Popen(
                [
                    fistar_path, unpacked_fname,
                    '--comment',
                    '--flux-threshold', '3000',
                    '--sort', 'flux',
                    '--format', 'id,x,y,s,d,k,flux,bg,s/n'
                ],
                stdout=PIPE
            ) as fistar:
                extracted = parse_result(fistar.stdout)
        print(f'Got:\n{extracted!r}')

        self.assertTrue(
            (
                expected.columns.size == extracted.columns.size
                and
                (expected.columns == extracted.columns).all()
            ),
            f'Source extraction of XO-1 image produced different columns than '
            f'expected:\n\t{extracted.columns!r}\n\tinstead of\n\t'
            f'{expected.columns!r}'
        )
        for column in expected.columns:
            expected_col = expected[column]
            extracted_col = extracted[column]
            self.assertTrue(
                expected_col.dtype == extracted_col.dtype,
                f'Column types mismatch: {extracted_col.dtype!r} instead of '
                f'{expected_col.dtype!r}'
            )
            mismatch_message = (
                f'Column {column!r} mismatch:\n\t{extracted_col!r}\n\t'
                f'instead of\n\t{expected_col!r}'
            )

            if expected_col.dtype.kind == 'i':
                self.assertTrue(expected_col.equals(extracted_col),
                                mismatch_message)
            else:
                assert expected_col.dtype.kind == 'f'
            for expected_val, extracted_val in zip(expected_col, extracted_col):
                self.assertApprox(expected_val, extracted_val, mismatch_message)


if __name__ == '__main__':
    unittest.main()
