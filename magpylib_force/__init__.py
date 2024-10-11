"""
The magpylib-force package sits on top of magpylib for force
computation between Magpylib source objects. API and examples
can be found in the Magylib documentation.
"""

from magpylib_force.force import getFT

# module level dunders
__version__ = "0.1.10"
__author__ = "SAL"
__all__ = [
    "getFT",
]
