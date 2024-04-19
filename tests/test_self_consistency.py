import numpy as np
import magpylib as magpy
from magpylib_force.force import getFTcube
from magpylib_force.force import getFTwire


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
#     cube1.mesh = (20,10,10)
#     cube2.mesh = (10,10,20)

#     F1,_ = getFTcube(cube1, cube2)
#     F2,_ = getFTcube(cube2, cube1)
#     errF = 2*(np.linalg.norm(F1+F2)/np.linalg.norm(F1-F2))
#     assert errF < 1e-4

#     ureg = cube1.dimension._REGISTRY
#     Tanch = np.array([0,0,0]) * ureg.meter
#     F1,T1 = getFTcube(cube1, cube2, Tanch=Tanch)
#     F2,T2 = getFTcube(cube2, cube1, Tanch=Tanch)

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
    cube1.mesh = (20,10,10)
    cube2.mesh = (10,10,20)

    #F1,_ = getFTcube(cube1, cube2)
    #F2,_ = getFTcube(cube2, cube1)
    #errF = 2*(np.linalg.norm(F1+F2)/np.linalg.norm(F1-F2))
    #assert errF < 1e-4

    F1,T1 = getFTcube(cube1, cube2, anchor=np.array((0,0,0)))
    F2,T2 = getFTcube(cube2, cube1, anchor=np.array((0,0,0)))

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
#     wire1.mesh=(200)
#     wire2 = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=([(-1,-1,.5), (-1,1,1), (1,1,2), (1,-1,1), (-1,-1,.5)], "m"),
#     )
#     wire2.mesh=(200)
#     F1a,_ = getFTwire(wire1, wire2)
#     F2a,_ = getFTwire(wire2, wire1)
#     errFa = np.linalg.norm(F1a+F2a) / np.linalg.norm(F1a-F2a)
#     assert errFa < 1e-5

#     ureg = wire1.current._REGISTRY
#     F1b,T1b = getFTwire(wire1, wire2, Tanch=np.array((0,0,0))*ureg.meter)
#     F2b,T2b = getFTwire(wire2, wire1, Tanch=np.array((0,0,0))*ureg.meter)
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
    wire1.mesh=(200)
    wire2 = magpy.current.Polyline(
        current=1,
        vertices=[(-1,-1,.5), (-1,1,1), (1,1,2), (1,-1,1), (-1,-1,.5)],
    )
    wire2.mesh=(200)
    F1a,_ = getFTwire(wire1, wire2)
    F2a,_ = getFTwire(wire2, wire1)
    errFa = np.linalg.norm(F1a+F2a) / np.linalg.norm(F1a-F2a)
    assert errFa < 1e-5

    F1b,T1b = getFTwire(wire1, wire2, anchor=np.array((0,0,0)))
    F2b,T2b = getFTwire(wire2, wire1, anchor=np.array((0,0,0)))
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
#     magnet.mesh=(30,20,10)
#     wire = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=([(-1,-1,.5), (-1,1,1), (1,1,.1), (1,-1,1), (-1,-1,.5)], "m"),
#     )
#     wire.mesh=(200)
#     F1a,_ = getFTcube(wire, magnet)
#     F2a,_ = getFTwire(magnet, wire)
#     errFa = np.linalg.norm(F1a+F2a) / np.linalg.norm(F1a-F2a)*2
#     assert errFa < 1e-5

#     ureg = wire.current._REGISTRY
#     F1b,T1b = getFTcube(wire, magnet, Tanch=np.array((0,0,0))*ureg.meter)
#     F2b,T2b = getFTwire(magnet, wire, Tanch=np.array((0,0,0))*ureg.meter)
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
    magnet.mesh=(30,20,10)
    wire = magpy.current.Polyline(
        current=1,
        vertices=[(-1,-1,.5), (-1,1,1), (1,1,.1), (1,-1,1), (-1,-1,.5)],
    )
    wire.mesh=(200)
    #F1a,_ = getFTcube(wire, magnet)
    #F2a,_ = getFTwire(magnet, wire)
    #errFa = np.linalg.norm(F1a+F2a) / np.linalg.norm(F1a-F2a)*2
    #assert errFa < 1e-5

    F1b,T1b = getFTcube(wire, magnet, anchor=np.array((0,0,0)))
    F2b,T2b = getFTwire(magnet, wire, anchor=np.array((0,0,0)))
    errFb = np.linalg.norm(F1b+F2b) / np.linalg.norm(F1b-F2b)*2
    assert errFb < 1e-5
    errTb = np.linalg.norm(T1b+T2b) / np.linalg.norm(T1b-T2b)*2
    assert errTb < 1e-5
