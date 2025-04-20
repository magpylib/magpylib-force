"""
Helper functions
"""

from __future__ import annotations

import warnings

import numpy as np
from magpylib._src.obj_classes.class_current_Circle import Circle
from magpylib._src.obj_classes.class_current_Polyline import Polyline
from magpylib._src.obj_classes.class_magnet_Cuboid import Cuboid
from magpylib._src.obj_classes.class_magnet_Cylinder import Cylinder
from magpylib._src.obj_classes.class_magnet_CylinderSegment import CylinderSegment
from magpylib._src.obj_classes.class_magnet_Sphere import Sphere


def check_input_anchor(anchor):
    """
    check and format anchor input
    """
    if anchor is None:
        warnings.warn(
            "No anchor was handed to getFT. This results in incorrect "
            "torque computation by excluding force contribution to torque.",
            stacklevel=2,
        )
        return None
    if isinstance(anchor, list | tuple):
        anchor = np.array(anchor)
    if not isinstance(anchor, np.ndarray):
        msg = "Anchor input must be list tuple or array."
        raise ValueError(msg)
    if anchor.shape != (3,):
        msg = "Anchor input must have shape (3,)."
        raise ValueError(msg)
    return anchor


def check_input_targets(targets):
    """check and format targets input"""
    if not isinstance(targets, list):
        targets = [targets]
    for t in targets:
        if not isinstance(
            t, Cuboid | Polyline | Sphere | Cylinder | CylinderSegment | Circle
        ):
            msg = (
                "Bad `targets` input for getFT."
                " `targets` can only be Cuboids, Polylines, Spheres, Cylinders, "
                " CylinderSegments, and Circles."
                f" Instead received type {type(t)} target."
            )
            raise ValueError(msg)
        if not hasattr(t, "meshing"):
            msg = (
                "Missing input for getFT `targets`."
                " `targets` must have the `meshing` parameter set."
            )
            raise ValueError(msg)
        if not isinstance(t, Polyline) and np.isscalar(t.meshing) and t.meshing < 20:
            warnings.warn(
                "Input parameter `meshing` is set to a low value which will result in "
                "inaccurate computation of force and torque.",
                stacklevel=2,
            )
    return targets
