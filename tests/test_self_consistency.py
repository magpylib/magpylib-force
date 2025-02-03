import warnings
import numpy as np
import magpylib as magpy
from magpylib_force.force import getFT


# def test_consistency_cube_cube():
#     """
#     consistency between 2 cuboids with aligned axes
#     """
#     cube1 = magpy.magnet.Cuboid(
#         dimension=((2,1,1), "mm"),
#         polarization=((1,2,3), "T"),
#     )
#     cube2 = magpy.magnet.Cuboid(
#         dimension=((1,1,2), "mm"),
#         polarization=((3,2,1), "T"),
#         position=((2,2,.5), "mm")
#     )
#     cube1.meshing = (20,10,10)
#     cube2.meshing = (10,10,20)

#     F1,_ = getFT(cube1, cube2)
#     F2,_ = getFT(cube2, cube1)
#     errF = 2*(np.linalg.norm(F1+F2)/np.linalg.norm(F1-F2))
#     assert errF < 1e-4

#     ureg = cube1.dimension._REGISTRY
#     Tanch = np.array([0,0,0]) * ureg.meter
#     F1,T1 = getFT(cube1, cube2, Tanch=Tanch)
#     F2,T2 = getFT(cube2, cube1, Tanch=Tanch)

#     errF = 2*(np.linalg.norm(F1+F2)/np.linalg.norm(F1-F2))
#     errT = 2*(np.linalg.norm(T1+T2)/np.linalg.norm(T1-T2))
#     assert errF < 1e-4
#     assert errT < 1e-4

def test_consistency_cube_cube():
    """
    consistency between 2 cuboids with aligned axes
    """
    cube1 = magpy.magnet.Cuboid(
        dimension=(2,1,1),
        polarization=(1,2,3),
    )
    cube2 = magpy.magnet.Cuboid(
        dimension=(1,1,2),
        polarization=(3,2,1),
        position=(2,2,.5)
    )
    cube1.meshing = (20,10,10)
    cube2.meshing = (10,10,20)

    #F1,_ = getFT(cube1, cube2)
    #F2,_ = getFT(cube2, cube1)
    #errF = 2*(np.linalg.norm(F1+F2)/np.linalg.norm(F1-F2))
    #assert errF < 1e-4

    F1,T1 = getFT(cube1, cube2, anchor=np.array((0,0,0)))
    F2,T2 = getFT(cube2, cube1, anchor=np.array((0,0,0)))

    errF = 2*(np.linalg.norm(F1+F2)/np.linalg.norm(F1-F2))
    errT = 2*(np.linalg.norm(T1+T2)/np.linalg.norm(T1-T2))
    assert errF < 1e-4
    assert errT < 1e-4


# def test_consistency_loop_loop():
#     """
#     consistency between 2 arbitrary current loops
#     """
#     wire1 = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=([(1,0,0),(0,1,.2), (-1,0,0), (0,-1,.5), (1,0,0)], "m"),
#     )
#     wire1.meshing=(200)
#     wire2 = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=([(-1,-1,.5), (-1,1,1), (1,1,2), (1,-1,1), (-1,-1,.5)], "m"),
#     )
#     wire2.meshing=(200)
#     F1a,_ = getFT(wire1, wire2)
#     F2a,_ = getFT(wire2, wire1)
#     errFa = np.linalg.norm(F1a+F2a) / np.linalg.norm(F1a-F2a)
#     assert errFa < 1e-5

#     ureg = wire1.current._REGISTRY
#     F1b,T1b = getFT(wire1, wire2, Tanch=np.array((0,0,0))*ureg.meter)
#     F2b,T2b = getFT(wire2, wire1, Tanch=np.array((0,0,0))*ureg.meter)
#     errFb = np.linalg.norm(F1b+F2b) / np.linalg.norm(F1b-F2b)
#     errTb = np.linalg.norm(T1b+T2b) / np.linalg.norm(T1b-T2b)
#     assert errFb < 1e-5
#     assert errTb < 1e-4

def test_consistency_loop_loop():
    """
    consistency between 2 arbitrary current loops
    """
    wire1 = magpy.current.Polyline(
        current=1,
        vertices=[(1,0,0),(0,1,.2), (-1,0,0), (0,-1,.5), (1,0,0)],
    )
    wire1.meshing=(200)
    wire2 = magpy.current.Polyline(
        current=1,
        vertices=[(-1,-1,.5), (-1,1,1), (1,1,2), (1,-1,1), (-1,-1,.5)],
    )
    wire2.meshing=(200)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        F1a,_ = getFT(wire1, wire2)
        F2a,_ = getFT(wire2, wire1)
    errFa = np.linalg.norm(F1a+F2a) / np.linalg.norm(F1a-F2a)
    assert errFa < 1e-5

    F1b,T1b = getFT(wire1, wire2, anchor=np.array((0,0,0)))
    F2b,T2b = getFT(wire2, wire1, anchor=np.array((0,0,0)))
    errFb = np.linalg.norm(F1b+F2b) / np.linalg.norm(F1b-F2b)
    errTb = np.linalg.norm(T1b+T2b) / np.linalg.norm(T1b-T2b)
    assert errFb < 1e-5
    assert errTb < 1e-4


# def test_consistency_cube_loop():
#     """
#     consistency between a current loop and a cuboid magnet
#     """
#     magnet = magpy.magnet.Cuboid(
#         polarization=((1,2,3), "T"),
#         dimension=((.6,.4,.2), "m"),
#     )
#     magnet.meshing=(30,20,10)
#     wire = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=([(-1,-1,.5), (-1,1,1), (1,1,.1), (1,-1,1), (-1,-1,.5)], "m"),
#     )
#     wire.meshing=(200)
#     F1a,_ = getFT(wire, magnet)
#     F2a,_ = getFT(magnet, wire)
#     errFa = np.linalg.norm(F1a+F2a) / np.linalg.norm(F1a-F2a)*2
#     assert errFa < 1e-5

#     ureg = wire.current._REGISTRY
#     F1b,T1b = getFT(wire, magnet, Tanch=np.array((0,0,0))*ureg.meter)
#     F2b,T2b = getFT(magnet, wire, Tanch=np.array((0,0,0))*ureg.meter)
#     errFb = np.linalg.norm(F1b+F2b) / np.linalg.norm(F1b-F2b)*2
#     assert errFb < 1e-5
#     errTb = np.linalg.norm(T1b+T2b) / np.linalg.norm(T1b-T2b)*2
#     assert errTb < 1e-5
def test_consistency_cube_loop():
    """
    consistency between a current loop and a cuboid magnet
    """
    magnet = magpy.magnet.Cuboid(
        polarization=(1,2,3),
        dimension=(.6,.4,.2),
    )
    magnet.meshing=(30,20,10)
    wire = magpy.current.Polyline(
        current=1,
        vertices=[(-1,-1,.5), (-1,1,1), (1,1,.1), (1,-1,1), (-1,-1,.5)],
    )
    wire.meshing=(200)
    #F1a,_ = getFT(wire, magnet)
    #F2a,_ = getFT(magnet, wire)
    #errFa = np.linalg.norm(F1a+F2a) / np.linalg.norm(F1a-F2a)*2
    #assert errFa < 1e-5

    F1b,T1b = getFT(wire, magnet, anchor=np.array((0,0,0)))
    F2b,T2b = getFT(magnet, wire, anchor=np.array((0,0,0)))
    errFb = np.linalg.norm(F1b+F2b) / np.linalg.norm(F1b-F2b)*2
    assert errFb < 1e-5
    errTb = np.linalg.norm(T1b+T2b) / np.linalg.norm(T1b-T2b)*2
    assert errTb < 1e-5


def test_consistency_cylinder_polyline():
    """
    check backward and forward for cylinder and polyline
    """
    loop = magpy.current.Polyline(
    vertices=((-3,0,0), (0,-3,0), (3,0,0), (0,3,0), (-3,0,0)),
    current=1000,
    position=(0,0,-1),
    )
    loop.rotate_from_angax(-45, (1,1,1))
    loop.meshing=100

    cyl = magpy.magnet.Cylinder(dimension=(2,1), polarization=(1,2,3))
    cyl.meshing=200

    ft1 = getFT(cyl, loop, anchor=cyl.position)
    ft2 = getFT(loop, cyl, anchor=cyl.position)

    assert np.max(abs((ft1+ft2)/(ft1-ft2))) < 0.001

    cyl.rotate_from_angax(33, (1,2,3))

    ft1 = getFT(cyl, loop, anchor=cyl.position)
    ft2 = getFT(loop, cyl, anchor=cyl.position)
    assert np.max(abs((ft1+ft2)/(ft1-ft2))) < 0.001


def test_consistency_cylinder_segment_cuboid():
    """
    check backward-forward with CylinderSegment and Cuboid
    """
    cube = magpy.magnet.Cuboid(
    dimension=(1,1,2),
    polarization=(-1,-2,-1)
    )
    cube.meshing=(8,8,4)

    dims = [
        (2,3,1,-20,120),
        (2,4,2,10,50),
        (3,4,3,50,100)
    ]
    pols = [
        (3,2,1),
        (1,2,3),
        (0,0,1),
    ]
    for dim,pol in zip(dims,pols):
        cyls = magpy.magnet.CylinderSegment(dimension=dim, polarization=pol)
        cyls.rotate_from_angax(-45, (1,1,1))
        cyls.meshing=300
        ft1 = getFT(cyls, cube, anchor=(0,0,0))
        ft2 = getFT(cube, cyls, anchor=(0,0,0))

        assert np.amax(abs((ft1+ft2)/(ft1-ft2))) < 0.09


def test_consistency_polyline_circle():
    """
    compare Polyline solution to circle solution
    """

    src = magpy.magnet.Sphere(diameter=1, polarization=(1,2,3), position=(0,0,-1))

    # circle
    loop1 = magpy.current.Circle(diameter=3, current=123)
    loop1.meshing=500
    # polyline
    rr = loop1.diameter/2
    ii = loop1.current
    phis = np.linspace(0,2*np.pi,500)
    verts = [(rr*np.cos(p), rr*np.sin(p), 0) for p in phis]
    loop2 = magpy.current.Polyline(current=ii, vertices=verts)
    loop2.meshing=1

    F1,T1 = getFT(src, loop1, anchor=(0,0,0))
    F2,T2 = getFT(src, loop2, anchor=(0,0,0))
    assert abs(np.linalg.norm(F1-F2)/np.linalg.norm(F1+F2)) < 1e-7
    assert abs(np.linalg.norm(T1-T2)/np.linalg.norm(T1+T2)) < 1e-7

    loop1.move((1.123,2.321,.123))
    loop2.move((1.123,2.321,.123))
    F1,T1 = getFT(src, loop1, anchor=(0,0,0))
    F2,T2 = getFT(src, loop2, anchor=(0,0,0))
    assert abs(np.linalg.norm(F1-F2)/np.linalg.norm(F1+F2)) < 1e-7
    assert abs(np.linalg.norm(T1-T2)/np.linalg.norm(T1+T2)) < 1e-7

    loop1.rotate_from_angax(20, 'x')
    loop2.rotate_from_angax(20, 'x')
    F1,T1 = getFT(src, loop1, anchor=(0,0,0))
    F2,T2 = getFT(src, loop2, anchor=(0,0,0))
    assert abs(np.linalg.norm(F1-F2)/np.linalg.norm(F1+F2)) < 1e-7
    assert abs(np.linalg.norm(T1-T2)/np.linalg.norm(T1+T2)) < 1e-7


def test_consistency_sphere_dipole():
    """
    force on sphere and dipole should be the same in nearly homogeneous field
    """

    src = magpy.current.Circle(diameter=10, current=123)
    pos = (0,0,0)
    diameter = 0.5
    magnetization_sphere = np.array((1e6,2e6,3e6))
    moment_dipole = magnetization_sphere*diameter**3/6*np.pi

    sphere = magpy.magnet.Sphere(diameter=diameter, magnetization=magnetization_sphere, position=pos)
    sphere.meshing = 100

    dipole = magpy.misc.Dipole(position=pos, moment=moment_dipole)

    FT_sphere = getFT(src, sphere, anchor=(0,0,0))
    FT_dipole = getFT(src, dipole, anchor=(0,0,0))

    np.testing.assert_allclose(FT_sphere, FT_dipole, rtol=1e-6, atol=1e-6)
