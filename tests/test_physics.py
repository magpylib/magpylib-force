import numpy as np
import magpylib as magpy
from magpylib_force.force import getFT


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
#     cloop.meshing = 1

#     # homogeneous field
#     def func(field, observers):
#         return np.zeros_like(observers, dtype=float) + np.array((1,0,0))
#     hom = magpy.misc.CustomSource(field_func=func)

#     # without anchor
#     F,T = getFT(hom, cloop, Tanch=None)
#     assert np.amax(F.magnitude) < 1e-14
#     assert T.magnitude == 0

#     # with anchor
#     F,T = getFT(hom, cloop, Tanch=cloop.position)
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
#     rloop.meshing = 10

#     # without anchor
#     F,T = getFT(hom, rloop, Tanch=None)
#     assert np.amax(F.magnitude) < 1e-14
#     assert T.magnitude == 0

#     # with anchor
#     F,T = getFT(hom, rloop, Tanch=rloop.position)
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
    cloop.meshing = 1

    # homogeneous field
    def func(field, observers):
        return np.zeros_like(observers, dtype=float) + np.array((1,0,0))
    hom = magpy.misc.CustomSource(field_func=func)

    # without anchor
    F,T = getFT(hom, cloop, anchor=None)
    assert np.amax(abs(F)) < 1e-14
    assert np.amax(abs(T)) == 0

    # with anchor
    F,T = getFT(hom, cloop, anchor=cloop.position)
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
    rloop.meshing = 10

    # without anchor
    F,T = getFT(hom, rloop, anchor=None)
    assert np.amax(abs(F)) < 1e-14
    assert np.amax(abs(T)) == 0

    # with anchor
    F,T = getFT(hom, rloop, anchor=rloop.position)
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
#     tgt.meshing = 1000

#     F,_ = getFT(src, tgt)

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
    wire2.meshing = 1000

    F,_ = getFT(wire1, wire2)

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
#     tgt.meshing = 1000

#     ureg=src.current._REGISTRY
#     F,T = getFT(src, tgt)

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
    wire2.meshing = 1000

    F,T = getFT(wire1, wire2)

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
    cube.meshing=(20,20,20)
    nn=150
    currs = []
    for z in np.linspace(-1,1,nn):
        loop = magpy.current.Polyline(
            current=1e7/4/np.pi*2/nn,
            vertices=((-1,-1,z),(1,-1,z),(1,1,z),(-1,1,z),(-1,-1,z)),
            position=(1,2,3)
        )
        loop.meshing = nn
        currs.append(loop)

    F1,T1 = getFT(gen, cube, anchor=np.array((0,0,0)))
    F2,T2 = np.sum(getFT(gen, currs, anchor=np.array((0,0,0))), axis=0)

    assert np.amax(abs((F1-F2)/(F1+F2)*2))<1e-2
    assert np.amax(abs((T1-T2)/(T1+T2)*2))<1e-2


def test_sphere_cube_at_distance():
    """
    A sphere and a cuboid with similar volume should see a similar torque and force
    at a distance
    """
    source = magpy.magnet.Sphere(diameter=1, polarization=(1,2,3))

    J = (3,2,1)
    pos = (5,-7,11)

    cube = magpy.magnet.Cuboid(
        dimension=(1,1,1),
        polarization=J,
        position=pos
    )
    cube.meshing = (2,2,2)

    sphere = magpy.magnet.Sphere(
        diameter=(6/np.pi)**(1/3),
        polarization=J,
        position=pos,
    )
    sphere.meshing=2

    FT = getFT(source, [cube, sphere], anchor=(0,0,0))

    errF = (FT[0,0]-FT[1,0])/np.linalg.norm(FT[0,0])
    errT = (FT[0,1]-FT[1,1])/np.linalg.norm(FT[0,1])

    assert max(abs(errF)) < 1e-5
    assert max(abs(errT)) < 1e-5
