"""Convenience for getting shared library filename for debugging purposes."""
from os import path

print(path.join(path.dirname(__file__), 'libastrowisp.dll'))
