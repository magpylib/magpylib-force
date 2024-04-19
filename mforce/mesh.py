import numpy as np
import itertools
import math as m
import pint

def mesh_cuboid(object):
    """
    splits cuboid into given mesh
    returns grid positions relative to cuboid position
    """
    a,b,c = object.dimension
    n1,n2,n3 = object.mesh
    xs = np.linspace(-a/2, a/2, n1+1)
    ys = np.linspace(-b/2, b/2, n2+1)
    zs = np.linspace(-c/2, c/2, n3+1)

    dx = xs[1] - xs[0] if len(xs)>1 else a
    dy = ys[1] - ys[0] if len(ys)>1 else b
    dz = zs[1] - zs[0] if len(zs)>1 else c

    xs_cent = xs[:-1] + dx/2 if len(xs)>1 else xs + dx/2
    ys_cent = ys[:-1] + dy/2 if len(ys)>1 else ys + dy/2
    zs_cent = zs[:-1] + dz/2 if len(zs)>1 else zs + dz/2

    if isinstance(a, pint.Quantity):
        permutas = np.array(list(itertools.product(xs_cent.magnitude, ys_cent.magnitude, zs_cent.magnitude)))
        return permutas * xs_cent.units

    return np.array(list(itertools.product(xs_cent, ys_cent, zs_cent)))




def mesh_cuboid_old2(object, verbose=False):
    """
    Splits up the object volume into a regular grid of small rectangular cells.
    The returned grid points lie in the center of these cells.

    Parameters
    ----------
    object: Magpylib source object with cuboid geometry
        Must have the parameter `dimension` with shape (3,).
        Must have the parameter `mesh` which is an int n or a triplet (n1,n2,n3).
        If mesh is int n, the cells are created aiming to achieve cubic aspect
        ratios. The resulting number of cells lies close to n. If input is
        triplet (n1,n2,n3), the mesh is created with those splittings.
    verbose: bool, default=False
        Print resulting mesh parameters

    Returns
    -------
    grid positions: np.ndarray shape (m,3)
    """
    a,b,c = object.dimension
    splitting = object.mesh
    if isinstance(splitting, (float, int)):
        x = (a*b*c/splitting)**(1/3)
        n1 = m.ceil(a/x)
        n2 = m.ceil(b/x)
        n3 = m.ceil(c/x)
    else:
        n1,n2,n3 = np.array(splitting)+1
    xs = np.linspace(-a/2, a/2, n1)
    ys = np.linspace(-b/2, b/2, n2)
    zs = np.linspace(-c/2, c/2, n3)

    dx = xs[1] - xs[0] if len(xs)>1 else a
    dy = ys[1] - ys[0] if len(ys)>1 else b
    dz = zs[1] - zs[0] if len(zs)>1 else c

    xs_cent = xs[:-1] + dx/2 if len(xs)>1 else xs + dx/2
    ys_cent = ys[:-1] + dy/2 if len(ys)>1 else ys + dy/2
    zs_cent = zs[:-1] + dz/2 if len(zs)>1 else zs + dz/2

    permutas = np.array(list(itertools.product(xs_cent, ys_cent, zs_cent)))

    if verbose:
        print('####### mesh statistics #######')
        print(f'- No. cells: {(n1-1)*(n2-1)*(n3-1)}')
        print(f'- Splitting: {n1-1} x {n2-1} x {n3-1}')
        print(f'- Cell dim : {np.round(dx, 4)} x {np.round(dy, 4)} x {np.round(dz, 4)}')

    offset = object.position
    rotation = object.orientation

    return rotation.apply(permutas) + offset




def cuboid_data(center, size):
    """
    Create a data array for cuboid plotting.
    ============= ================================================
    Argument      Description
    ============= ================================================
    center        center of the cuboid, triple
    size          size of the cuboid, triple, (x_length,y_width,z_height)
    :type size: tuple, numpy.array, list
    :param size: size of the cuboid, triple, (x_length,y_width,z_height)
    :type center: tuple, numpy.array, list
    :param center: center of the cuboid, triple, (x,y,z)
    """

    # suppose axis direction: x: to left; y: to inside; z: to upper
    # get the (left, outside, bottom) point
    o = [a - b / 2 for a, b in zip(center, size)]
    # get the length, width, and height
    l, w, h = size
    x = [[o[0], o[0] + l, o[0] + l, o[0], o[0]],  # x coordinate of points in bottom surface
         [o[0], o[0] + l, o[0] + l, o[0], o[0]],  # x coordinate of points in upper surface
         [o[0], o[0] + l, o[0] + l, o[0], o[0]],  # x coordinate of points in outside surface
         [o[0], o[0] + l, o[0] + l, o[0], o[0]]]  # x coordinate of points in inside surface
    y = [[o[1], o[1], o[1] + w, o[1] + w, o[1]],  # y coordinate of points in bottom surface
         [o[1], o[1], o[1] + w, o[1] + w, o[1]],  # y coordinate of points in upper surface
         [o[1], o[1], o[1], o[1], o[1]],          # y coordinate of points in outside surface
         [o[1] + w, o[1] + w, o[1] + w, o[1] + w, o[1] + w]]    # y coordinate of points in inside surface
    z = [[o[2], o[2], o[2], o[2], o[2]],                        # z coordinate of points in bottom surface
         [o[2] + h, o[2] + h, o[2] + h, o[2] + h, o[2] + h],    # z coordinate of points in upper surface
         [o[2], o[2], o[2] + h, o[2] + h, o[2]],                # z coordinate of points in outside surface
         [o[2], o[2], o[2] + h, o[2] + h, o[2]]]                # z coordinate of points in inside surface
    return x, y, z
