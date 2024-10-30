import warnings
import numpy as np
import magpylib as magpy
from magpylib_force.force import getFT
from scipy.special import ellipe
from scipy.special import ellipk


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
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        F,T = getFT(hom, cloop, anchor=None)
    assert np.amax(abs(F)) < 1e-14
    assert np.amax(abs(T)) == 0

    # with anchor
    F,T = getFT(hom, cloop, anchor=cloop.position)
    T*=-1 #bad sign at initial test design
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
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        F,T = getFT(hom, rloop, anchor=None)
    assert np.amax(abs(F)) < 1e-14
    assert np.amax(abs(T)) == 0

    # with anchor
    F,T = getFT(hom, rloop, anchor=rloop.position)
    T*=-1 #bad sign at initial test design
    assert np.amax(abs(F)) < 1e-14
    assert abs(T[0]) < 1e-14
    assert abs(T[1] + 4 ) < 1e-3
    assert abs(T[2]) < 1e-14


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

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        F,_ = getFT(wire1, wire2)

    Fanalytic = 2*magpy.mu_0/4/np.pi*2000

    assert abs(F[0]) < 1e-14
    assert abs(F[1]) < 1e-14
    assert abs((F[2]+Fanalytic)/Fanalytic) < 1e-3


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

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        F,_ = getFT(wire1, wire2)

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


def test_physics_at_distance():
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
    cube.meshing = 100

    sphere = magpy.magnet.Sphere(
        diameter=(6/np.pi)**(1/3),
        polarization=J,
        position=pos,
    )
    sphere.meshing=100

    cyl = magpy.magnet.Cylinder(
        dimension=(2*np.sqrt(1/np.pi),1),
        polarization=J,
        position=pos,
    )
    cyl.meshing=100

    FT = getFT(source, [cube, sphere, cyl], anchor=(0,0,0))
    print(FT.shape)

    for i in range(1,3):
        errF = abs(np.linalg.norm(FT[0,0]-FT[i,0]) / np.linalg.norm(FT[0,0]+FT[i,0]))
        errT = abs(np.linalg.norm(FT[0,1]-FT[i,1]) / np.linalg.norm(FT[0,1]+FT[i,1]))
        assert errF<1e-4
        assert errF<1e-4


def test_physics_torque_sign():
    """ make sure that torque sign is in the right direction"""

    # Cuboid -> Cuboid
    mag1 = magpy.magnet.Cuboid(position=(2,0,0), polarization=(1,0,0), dimension=(2,1,1))
    mag2 = magpy.magnet.Cuboid(position=(-2,0,0), polarization=(1,0,0), dimension=(2,1,1))

    mag1.rotate_from_angax(15, "y")
    mag1.meshing=(3,3,3)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _,T = getFT(mag2, mag1)

    assert T[1] < 0

    # Cuboid -> Polyline
    mag = magpy.magnet.Cuboid(polarization=(0,0,1), dimension=(1,1,2))

    ts = np.linspace(0,2*np.pi,10)
    verts = [(2*np.cos(t), 2*np.sin(t),0) for t in ts]
    loop = magpy.current.Polyline(vertices=verts, current=1)
    loop.rotate_from_angax(15, "y")

    loop.meshing=2

    _,T = getFT(mag, loop, anchor=(0,0,0))

    assert T[1] < 0


def test_physics_force_between_cocentric_loops():
    """
    compare the numerical solution against the analytical solution of the force between two
    cocentric current loops.
    See e.g. IEEE TRANSACTIONS ON MAGNETICS, VOL. 49, NO. 8, AUGUST 2013
    """
    z1, z2 = 0.123, 1.321
    i1, i2 = 3.2, 5.1
    r1, r2 = 1.2, 2.3

    # numerical solution
    loop1 = magpy.current.Circle(diameter=2*r1, current=i1, position=(0,0,z1))
    loop2 = magpy.current.Circle(diameter=2*r2, current=i2, position=(0,0,z2))
    loop2.meshing=1000
    F_num = getFT(loop1, loop2, anchor=(0,0,0))[0,2]

    # analytical solution
    k2 = 4*r1*r2 / ((r1+r2)**2+(z1-z2)**2)
    k = np.sqrt(k2)
    pf = magpy.mu_0*i1*i2*(z1-z2)*k / 4 / np.sqrt(r1*r2)
    F_ana = pf*( (2-k2)/(1-k2)*ellipe(k**2) - 2*ellipk(k**2) )

    assert abs((F_num - F_ana)/(F_num + F_ana)) < 1e-5
