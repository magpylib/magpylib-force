import numpy as np
import magpylib as magpy
from magpylib_force.force import getFTcube
from magpylib_force.force import getFTwire


# def test_physics_loop_torque():
#     """
#     for a current loop in a homogeneous field the following holds
#     F = 0
#     T = current * loop_surface * field_normal_component
#     """
#     # circular loop
#     ts = np.linspace(0,2*np.pi,300)
#     verts = [(np.sin(t), np.cos(t), 0) for t in ts]
#     cloop = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=(verts, "m"),
#     )
#     cloop.mesh = 1

#     # homogeneous field
#     def func(field, observers):
#         return np.zeros_like(observers, dtype=float) + np.array((1,0,0))
#     hom = magpy.misc.CustomSource(field_func=func)

#     # without anchor
#     F,T = getFTwire(hom, cloop, Tanch=None)
#     assert np.amax(F.magnitude) < 1e-14
#     assert T.magnitude == 0

#     # with anchor
#     F,T = getFTwire(hom, cloop, Tanch=cloop.position)
#     assert np.amax(abs(F.magnitude)) < 1e-14
#     assert abs(T[0].magnitude) < 1e-14
#     assert abs(T[1].magnitude - np.pi ) < 1e-3
#     assert abs(T[2].magnitude) < 1e-14

#     ##############################################################

#     # rectangular loop
#     verts = [(-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (-1,-1,0)]
#     rloop = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=(verts, "m"),
#     )
#     rloop.mesh = 10

#     # without anchor
#     F,T = getFTwire(hom, rloop, Tanch=None)
#     assert np.amax(F.magnitude) < 1e-14
#     assert T.magnitude == 0

#     # with anchor
#     F,T = getFTwire(hom, rloop, Tanch=rloop.position)
#     assert np.amax(abs(F.magnitude)) < 1e-14
#     assert abs(T[0].magnitude) < 1e-14
#     assert abs(T[1].magnitude + 4 ) < 1e-3
#     assert abs(T[2].magnitude) < 1e-14

def test_physics_loop_torque():
    """
    for a current loop in a homogeneous field the following holds
    F = 0
    T = current * loop_surface * field_normal_component
    """
    # circular loop
    ts = np.linspace(0,2*np.pi,300)
    verts = [(np.sin(t), np.cos(t), 0) for t in ts]
    cloop = magpy.current.Polyline(
        current=1,
        vertices=verts,
    )
    cloop.mesh = 1

    # homogeneous field
    def func(field, observers):
        return np.zeros_like(observers, dtype=float) + np.array((1,0,0))
    hom = magpy.misc.CustomSource(field_func=func)

    # without anchor
    F,T = getFTwire(hom, cloop, anchor=None)
    assert np.amax(F) < 1e-14
    assert T == 0

    # with anchor
    F,T = getFTwire(hom, cloop, anchor=cloop.position)
    assert np.amax(abs(F)) < 1e-14
    assert abs(T[0]) < 1e-14
    assert abs(T[1] - np.pi ) < 1e-3
    assert abs(T[2]) < 1e-14

    ##############################################################

    # rectangular loop
    verts = [(-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0), (-1,-1,0)]
    rloop = magpy.current.Polyline(
        current=1,
        vertices=verts,
    )
    rloop.mesh = 10

    # without anchor
    F,T = getFTwire(hom, rloop, anchor=None)
    assert np.amax(F) < 1e-14
    assert T == 0

    # with anchor
    F,T = getFTwire(hom, rloop, anchor=rloop.position)
    assert np.amax(abs(F)) < 1e-14
    assert abs(T[0]) < 1e-14
    assert abs(T[1] + 4 ) < 1e-3
    assert abs(T[2]) < 1e-14


# def test_physics_parallel_wires():
#     """
#     The force between straight infinite parallel wires is
#     F = 2*mu0/4/pi * i1*i2/r
#     """
#     src = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=([(-1000,0,0),(1000,0,0)], "m"),
#     )
#     tgt = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=([(-1000,0,0),(0,0,0),(1000,0,0)], "m"),
#         position=((0,0,1), "m"),
#     )
#     tgt.mesh = 1000

#     F,_ = getFTwire(src, tgt)

#     Fanalytic = 2*magpy.mu_0/4/np.pi*2000

#     assert abs(F[0].magnitude) < 1e-14
#     assert abs(F[1].magnitude) < 1e-14
#     assert abs((F[2].magnitude + Fanalytic)/Fanalytic) < 1e-3

def test_physics_parallel_wires():
    """
    The force between straight infinite parallel wires is
    F = 2*mu0/4/pi * i1*i2/r
    """
    wire1 = magpy.current.Polyline(
        current=1,
        vertices=[(-1000,0,0),(1000,0,0)],
    )
    wire2 = wire1.copy(position=(0,0,1))
    wire2.mesh = 1000

    F,_ = getFTwire(wire1, wire2)

    Fanalytic = 2*magpy.mu_0/4/np.pi*2000

    assert abs(F[0]) < 1e-14
    assert abs(F[1]) < 1e-14
    assert abs((F[2]+Fanalytic)/Fanalytic) < 1e-3


# def test_physics_perpendicular_wires():
#     """
#     The force between straight infinite perpendicular wires is 0
#     """
#     src = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=([(-1000,0,0),(1000,0,0)], "m"),
#     )
#     tgt = magpy.current.Polyline(
#         current=(1, "A"),
#         vertices=([(0,-1000,0),(0,0,0),(0,1000,0)], "m"),
#         position=((0,0,1), "m"),
#     )
#     tgt.mesh = 1000

#     ureg=src.current._REGISTRY
#     F,T = getFTwire(src, tgt)

#     assert np.max(abs(F.magnitude)) < 1e-14



def test_physics_perpendicular_wires():
    """
    The force between straight infinite perpendicular wires is 0
    """
    wire1 = magpy.current.Polyline(
        current=1,
        vertices=[(-1000,0,0),(1000,0,0)],
    )
    wire2 = magpy.current.Polyline(
        current=1,
        vertices=[(0,-1000,0),(0,0,0),(0,1000,0)],
        position=(0,0,1),
    )
    wire2.mesh = 1000

    F,T = getFTwire(wire1, wire2)

    assert np.max(abs(F)) < 1e-14


def test_cube_loop_replacement():
    """
    A magnet sees the same force as a respective current replacement
    """
    gen = magpy.magnet.Cuboid(
        polarization=(1,2,3),
        dimension=(1,1,1)
    )
    cube = magpy.magnet.Cuboid(
        polarization=(0,0,1),
        dimension=(2,2,2),
        position=(1,2,3)
    )
    cube.mesh=(20,20,20)
    nn=150
    currs = []
    for z in np.linspace(-1,1,nn):
        loop = magpy.current.Polyline(
            current=1e7/4/np.pi*2/nn,
            vertices=((-1,-1,z),(1,-1,z),(1,1,z),(-1,1,z),(-1,-1,z)),
            position=(1,2,3)
        )
        loop.mesh = nn
        currs.append(loop)

    F1,T1 = getFTcube(gen, cube, anchor=np.array((0,0,0)))
    F2,T2 = getFTwire(gen, currs, anchor=np.array((0,0,0)))
    F2 = np.sum(F2,axis=0)
    T2 = np.sum(T2,axis=0)

    assert np.amax(abs((F1-F2)/(F1+F2)*2))<1e-2
    assert np.amax(abs((T1-T2)/(T1+T2)*2))<1e-2
