"""
Copyright (c) 2025 Michael Ortner. All rights reserved.

magpylib-force: Python package extending the Magpylib library that enables force and torque computationss.
"""

from __future__ import annotations

from ._version import version as __version__

from magpylib_force.force import getFT

__all__ = ["__version__", "getFT"]

