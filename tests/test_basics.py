from __future__ import annotations

import magpylib as magpy
import numpy as np
from scipy.spatial.transform import Rotation as R

from magpylib_force import getFT


def test_rotation1():
    """
    test if rotated magnets give the correct result
    """
    s1 = magpy.magnet.Sphere(diameter=1, polarization=(1, 2, 3))

    c1 = magpy.magnet.Cuboid(
        dimension=(1, 1, 1),
        polarization=(0, 0, 1),
        position=(1, 2, 3),
    )
    c1.meshing = (5, 5, 5)  # regular rectangular mesh in Cuboid

    c2 = magpy.magnet.Cuboid(
        dimension=(1, 1, 1),
        polarization=(1, 0, 0),
        position=(1, 2, 3),
    )
    c2.meshing = (5, 5, 5)  # regular rectangular mesh in Cuboid
    c2.rotate_from_angax(-90, "y")

    FT = getFT(s1, [c1, c2], anchor=(0, 0, 0))

    np.testing.assert_allclose(FT[0], FT[1])


def test_rotation2():
    """
    test if rotated currents give the same result
    """
    s1 = magpy.magnet.Sphere(diameter=1, polarization=(1, 2, 3), position=(0, 0, -1))

    verts1 = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 0)]
    verts2 = [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1), (0, 0, 0)]

    c1 = magpy.current.Polyline(
        vertices=verts1,
        current=1,
    )
    c1.meshing = 15

    c2 = magpy.current.Polyline(
        vertices=verts2,
        current=1,
    )
    c2.meshing = 15
    c2.rotate_from_angax(-90, "x")

    FT = getFT(s1, [c1, c2], anchor=(0, 0, 0))

    np.testing.assert_allclose(FT[0], FT[1])


def test_orientation():
    """
    test if dipole with orientation gives same result as rotated magnetic moment
    """
    mm, md = np.array((0.976, 4.304, 2.055)), np.array((0.878, -1.527, 2.918))
    pm, pd = np.array((-1.248, 7.835, 9.273)), np.array((-2.331, 5.835, 0.578))

    magnet = magpy.magnet.Cuboid(position=pm, dimension=(1, 2, 3), polarization=mm)
    r = R.from_euler("xyz", (25, 65, 150), degrees=True)
    dipole1 = magpy.misc.Dipole(position=pd, moment=md, orientation=r)
    dipole2 = magpy.misc.Dipole(position=pd, moment=r.apply(md))

    F = getFT(magnet, [dipole1, dipole2], anchor=(0, 0, 0))

    np.testing.assert_allclose(F[0], F[1])
