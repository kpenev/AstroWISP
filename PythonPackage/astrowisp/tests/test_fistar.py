#!/usr/bin/env python3

"""Test source extraction using built-in fistar executable."""

from os import path
from subprocess import Popen, PIPE
from functools import partial
import platform

import unittest
from pandas import read_csv

from astrowisp import fistar_path
from astrowisp.utils.file_utilities import get_unpacked_fits
from astrowisp.tests.utilities import FloatTestCase

_test_data_dir = path.join(path.dirname(path.abspath(__file__)), "test_data")


class TestFistar(FloatTestCase):
    """Test cases for the fistar executable."""

    def test_xo1(self):
        """Check if extracting sources from XO-1 image matches expected."""

        parse_result = partial(read_csv, sep=r"\s+", comment="#", header=None)

        expected = parse_result(
            path.join(
                _test_data_dir,
                (
                    "XO1_test_img_applechip.fistar"
                    if platform.machine() == "arm64"
                    else "XO1_test_img.fistar"
                ),
            )
        )
        with get_unpacked_fits(
            path.join(_test_data_dir, "XO1_test_img.fits")
        ) as unpacked_fname:
            with Popen(
                [
                    fistar_path,
                    unpacked_fname,
                    "--comment",
                    "--flux-threshold",
                    "3000",
                    "--sort",
                    "flux",
                    "--format",
                    "id,x,y,s,d,k,flux,bg,s/n",
                ],
                stdout=PIPE,
            ) as fistar:
                extracted = parse_result(fistar.stdout)

        self.assertApproxPandas(
            expected, extracted, "Source extraction of XO-1 image"
        )


if __name__ == "__main__":
    unittest.main()
