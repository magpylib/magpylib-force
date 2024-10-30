"""
The magpylib-force package sits on top of magpylib for force
computation between Magpylib source objects. API and examples
can be found in the Magylib documentation.
"""

from magpylib_force.force import getFT

# module level dunders
__version__ = "0.3.1"
__author__ = "The Magpylib Team"
__all__ = [
    "getFT",
]
