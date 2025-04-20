from __future__ import annotations

import importlib.metadata

import magpylib_force as m


def test_version():
    assert importlib.metadata.version("magpylib_force") == m.__version__
