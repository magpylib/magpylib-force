import numpy as np
import itertools
import math as m
from magpylib._src.obj_classes.class_magnet_Cuboid import Cuboid
from magpylib._src.obj_classes.class_magnet_Sphere import Sphere
from magpylib._src.obj_classes.class_magnet_Cylinder import Cylinder
from magpylib._src.obj_classes.class_magnet_CylinderSegment import CylinderSegment

def mesh_target(object):
    """
    create mesh for target objects
    """
    if isinstance(object, Cuboid):
        return mesh_cuboid(object)
    elif isinstance(object, Sphere):
        return mesh_sphere(object)
    elif isinstance(object, Cylinder):
        return mesh_cylinder(object)
    raise RuntimeError("fktn `mesh_target`: should not be here!")


def mesh_sphere(object):
    """
    create sphere mesh from object meshing parameter
    """
    n = object.meshing
    dia = object.diameter
    a = -dia/2+dia/(2*n)
    b =  dia/2-dia/(2*n)
    c = n*1j
    mesh = np.mgrid[a:b:c, a:b:c, a:b:c].T.reshape(n**3,3)
    return mesh[np.linalg.norm(mesh, axis=1)<dia/2]


def mesh_cylinder(object):
    """
    create cylinder mesh from object meshing parameter
    """
    n = object.meshing
    dia, height = object.dimension

    a1 = -dia/2+dia/(2*n)
    b1 =  dia/2-dia/(2*n)
    a2 = -height/2+height/(2*n)
    b2 =  height/2-height/(2*n)

    if dia > height:
        c1 = n
        c2 = int(n/dia*height)
        if c2<3:
            print("Warning: bad cylinder mesh ratio. Increase meshing parameter.")
            c2 = 3
    else:
        c2 = n
        c1 = int(n/height*dia)
        if c1<3:
            print("Warning: bad cylinder mesh ratio. Increase meshing parameter.")
            c1 = 3

    mesh = np.mgrid[a1:b1:c1*1j, a1:b1:c1*1j, a2:b2:c2*1j].T.reshape(c1*c1*c2,3)
    mask = np.linalg.norm(mesh[:,:2], axis=1) < dia/2
    return mesh[mask]


def mesh_cylinder_segment(object):
    """
    create cylinder mesh from object meshing parameter
    """
    n = object.meshing
    r1, r2, h, phi1, phi2 = object.dimension

    a1 = -r2+r2/n
    b1 =  r2-r2/n
    a2 = -h/2 + h/(2*n)
    b2 =  h/2 - h/(2*n)

    dia = 2*r2
    if dia > h:
        c1 = n
        c2 = int(n/dia*h)
        c2 = 3 if c2<3 else c2
    else:
        c2 = n
        c1 = int(n/h*dia)
        c1 = 3 if c1<3 else c1

    mesh = np.mgrid[a1:b1:c1*1j, a1:b1:c1*1j, a2:b2:c2*1j].T.reshape(c1*c1*c2,3)
    import magpylib as magpy

    mask1 = np.linalg.norm(mesh[:,:2], axis=1) > r1
    mask2 = np.linalg.norm(mesh[:,:2], axis=1) < r2
    mask3 = np.arctan2(mesh[:,0], mesh[:,1]) > phi1/180*np.pi
    mask4 = np.arctan2(mesh[:,0], mesh[:,1]) < phi2/180*np.pi
    return mesh[mask1*mask2*mask3*mask4]


def mesh_cuboid(object):
    """
    splits cuboid into given mesh
    returns grid positions relative to cuboid position
    """
    a,b,c = object.dimension
    n1,n2,n3 = object.meshing
    xs = np.linspace(-a/2, a/2, n1+1)
    ys = np.linspace(-b/2, b/2, n2+1)
    zs = np.linspace(-c/2, c/2, n3+1)

    dx = xs[1] - xs[0] if len(xs)>1 else a
    dy = ys[1] - ys[0] if len(ys)>1 else b
    dz = zs[1] - zs[0] if len(zs)>1 else c

    xs_cent = xs[:-1] + dx/2 if len(xs)>1 else xs + dx/2
    ys_cent = ys[:-1] + dy/2 if len(ys)>1 else ys + dy/2
    zs_cent = zs[:-1] + dz/2 if len(zs)>1 else zs + dz/2

    # if isinstance(a, pint.Quantity):
    #     permutas = np.array(list(itertools.product(xs_cent.magnitude, ys_cent.magnitude, zs_cent.magnitude)))
    #     return permutas * xs_cent.units

    return np.array(list(itertools.product(xs_cent, ys_cent, zs_cent)))


# def mesh_cuboid_old2(object, verbose=False):
#     """
#     Splits up the object volume into a regular grid of small rectangular cells.
#     The returned grid points lie in the center of these cells.

#     Parameters
#     ----------
#     object: Magpylib source object with cuboid geometry
#         Must have the parameter `dimension` with shape (3,).
#         Must have the parameter `mesh` which is an int n or a triplet (n1,n2,n3).
#         If mesh is int n, the cells are created aiming to achieve cubic aspect
#         ratios. The resulting number of cells lies close to n. If input is
#         triplet (n1,n2,n3), the mesh is created with those splittings.
#     verbose: bool, default=False
#         Print resulting mesh parameters

#     Returns
#     -------
#     grid positions: np.ndarray shape (m,3)
#     """
#     a,b,c = object.dimension
#     splitting = object.meshing
#     if isinstance(splitting, (float, int)):
#         x = (a*b*c/splitting)**(1/3)
#         n1 = m.ceil(a/x)
#         n2 = m.ceil(b/x)
#         n3 = m.ceil(c/x)
#     else:
#         n1,n2,n3 = np.array(splitting)+1
#     xs = np.linspace(-a/2, a/2, n1)
#     ys = np.linspace(-b/2, b/2, n2)
#     zs = np.linspace(-c/2, c/2, n3)

#     dx = xs[1] - xs[0] if len(xs)>1 else a
#     dy = ys[1] - ys[0] if len(ys)>1 else b
#     dz = zs[1] - zs[0] if len(zs)>1 else c

#     xs_cent = xs[:-1] + dx/2 if len(xs)>1 else xs + dx/2
#     ys_cent = ys[:-1] + dy/2 if len(ys)>1 else ys + dy/2
#     zs_cent = zs[:-1] + dz/2 if len(zs)>1 else zs + dz/2

#     permutas = np.array(list(itertools.product(xs_cent, ys_cent, zs_cent)))

#     if verbose:
#         print('####### mesh statistics #######')
#         print(f'- No. cells: {(n1-1)*(n2-1)*(n3-1)}')
#         print(f'- Splitting: {n1-1} x {n2-1} x {n3-1}')
#         print(f'- Cell dim : {np.round(dx, 4)} x {np.round(dy, 4)} x {np.round(dz, 4)}')

#     offset = object.position
#     rotation = object.orientation

#     return rotation.apply(permutas) + offset




# def cuboid_data(center, size):
#     """
#     Create a data array for cuboid plotting.
#     ============= ================================================
#     Argument      Description
#     ============= ================================================
#     center        center of the cuboid, triple
#     size          size of the cuboid, triple, (x_length,y_width,z_height)
#     :type size: tuple, numpy.array, list
#     :param size: size of the cuboid, triple, (x_length,y_width,z_height)
#     :type center: tuple, numpy.array, list
#     :param center: center of the cuboid, triple, (x,y,z)
#     """

#     # suppose axis direction: x: to left; y: to inside; z: to upper
#     # get the (left, outside, bottom) point
#     o = [a - b / 2 for a, b in zip(center, size)]
#     # get the length, width, and height
#     l, w, h = size
#     x = [[o[0], o[0] + l, o[0] + l, o[0], o[0]],  # x coordinate of points in bottom surface
#          [o[0], o[0] + l, o[0] + l, o[0], o[0]],  # x coordinate of points in upper surface
#          [o[0], o[0] + l, o[0] + l, o[0], o[0]],  # x coordinate of points in outside surface
#          [o[0], o[0] + l, o[0] + l, o[0], o[0]]]  # x coordinate of points in inside surface
#     y = [[o[1], o[1], o[1] + w, o[1] + w, o[1]],  # y coordinate of points in bottom surface
#          [o[1], o[1], o[1] + w, o[1] + w, o[1]],  # y coordinate of points in upper surface
#          [o[1], o[1], o[1], o[1], o[1]],          # y coordinate of points in outside surface
#          [o[1] + w, o[1] + w, o[1] + w, o[1] + w, o[1] + w]]    # y coordinate of points in inside surface
#     z = [[o[2], o[2], o[2], o[2], o[2]],                        # z coordinate of points in bottom surface
#          [o[2] + h, o[2] + h, o[2] + h, o[2] + h, o[2] + h],    # z coordinate of points in upper surface
#          [o[2], o[2], o[2] + h, o[2] + h, o[2]],                # z coordinate of points in outside surface
#          [o[2], o[2], o[2] + h, o[2] + h, o[2]]]                # z coordinate of points in inside surface
#     return x, y, z
