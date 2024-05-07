import numpy as np
import magpylib as magpy
from magpylib_force import getFT


def test_rotation1():
    """
    test if rotated magnets give the correct result
    """
    s1 = magpy.magnet.Sphere(diameter=1, polarization=(1,2,3))

    c1 = magpy.magnet.Cuboid(
        dimension=(1,1,1),
        polarization=(0,0,1),
        position=(1,2,3),
    )
    c1.meshing = (5,5,5) # regular rectangular mesh in Cuboid

    c2 = magpy.magnet.Cuboid(
        dimension=(1,1,1),
        polarization=(1,0,0),
        position=(1,2,3),
    )
    c2.meshing = (5,5,5) # regular rectangular mesh in Cuboid
    c2.rotate_from_angax(-90, 'y')

    FT = getFT(s1, [c1, c2], anchor=(0,0,0))

    np.testing.assert_allclose(FT[0], FT[1])


def test_rotation2():
    """
    test if rotated currents give the same result
    """
    s1 = magpy.magnet.Sphere(diameter=1, polarization=(1,2,3), position=(0,0,-1))

    verts1 = [(0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,0)]
    verts2 = [(0,0,0), (1,0,0), (1,0,1), (0,0,1), (0,0,0)]

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

    FT = getFT(s1, [c1, c2], anchor=(0,0,0))

    np.testing.assert_allclose(FT[0], FT[1])
