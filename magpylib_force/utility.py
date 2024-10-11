"""
Helper functions
"""

import warnings
import numpy as np
from magpylib._src.obj_classes.class_magnet_Cuboid import Cuboid
from magpylib._src.obj_classes.class_current_Polyline import Polyline
from magpylib._src.obj_classes.class_magnet_Sphere import Sphere
from magpylib._src.obj_classes.class_magnet_Cylinder import Cylinder
from magpylib._src.obj_classes.class_magnet_CylinderSegment import CylinderSegment

def check_input_anchor(anchor):
    """
    check and format anchor input
    """
    if anchor is None:
        warnings.warn(
            "No anchor was handed to getFT. This results in incorrect "
            "torque computation by excluding force contribution to torque."
            )
        return None
    if isinstance(anchor, (list, tuple)):
        anchor = np.array(anchor)
    if not isinstance(anchor, np.ndarray):
        raise ValueError("Anchor input must be list tuple or array.")
    if anchor.shape != (3,):
        raise ValueError("Anchor input must have shape (3,).")
    return anchor


def check_input_targets(targets):
    """ check and format targets input """
    if not isinstance(targets, list):
        targets = [targets]
    for t in targets:
        if not isinstance(t, (Cuboid, Polyline, Sphere, Cylinder, CylinderSegment)):
            raise ValueError(
                "Bad `targets` input for getFT."
                " `targets` can only be Cuboids, Polylines, Spheres, Cylinders, and"
                " CylinderSegments."
                f" Instead receivd type {type(t)} target."
            )
        if not hasattr(t, "meshing"):
            raise ValueError(
                "Bad `targets` input for getFT."
                " `targets` must have the `meshing` parameter set."
            )
    return targets
