#!/usr/bin/env python3

"""List all modules imported by the python package."""

from sys import argv
import os
import fnmatch
from os import path
import re

if __name__ == "__main__":
    import_rex = re.compile(r"import|from \w* import")
    if len(argv) > 1:
        search_dir = argv[1]
    else:
        search_dir = path.join(
            path.dirname(path.abspath(__file__)), "PythonPackage"
        )
    for dirpath, dirnames, filenames in os.walk(search_dir):
        for fname in fnmatch.filter(filenames, "*.py"):
            fname = path.join(dirpath, fname)
            with open(fname, "r", encoding="ascii") as contents:
                line = contents.readline()
                while line:
                    if import_rex.match(line.strip()):
                        while line.strip().endswith("\\"):
                            line = (
                                line.rstrip()[:-1] + " " + contents.readline()
                            )
                        print(line.strip())
                        if line.strip().endswith("("):
                            while not line.strip().endswith(")"):
                                line = contents.readline()
                                print(line.strip())
                    line = contents.readline()
