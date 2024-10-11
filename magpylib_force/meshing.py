import numpy as np
import itertools
import math as m
from magpylib._src.obj_classes.class_magnet_Cuboid import Cuboid
from magpylib._src.obj_classes.class_magnet_Sphere import Sphere
from magpylib._src.obj_classes.class_magnet_Cylinder import Cylinder
from magpylib._src.obj_classes.class_magnet_CylinderSegment import CylinderSegment
from itertools import product


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
    elif isinstance(object, CylinderSegment):
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


def apportion_triple(triple, min_val=1, max_iter=30):
    """Apportion values of a triple, so that the minimum value `min_val` is respected
    and the product of all values remains the same.
    Example: apportion_triple([1,2,50], min_val=3)
    -> [ 2.99999999  3.         11.11111113]
    """
    triple = np.abs(np.array(triple, dtype=float))
    count = 0
    while any(n < min_val for n in triple) and count < max_iter:
        count += 1
        amin, amax = triple.argmin(), triple.argmax()
        factor = min_val / triple[amin]
        if triple[amax] >= factor * min_val:
            triple /= factor**0.5
            triple[amin] *= factor**1.5
    return triple

def cells_from_dimension(
    dim,
    target_elems,
    min_val=1,
    strict_max=False,
    parity=None,
):
    """Divide a dimension triple with a target scalar of elements, while apportioning
    the number of elements based on the dimension proportions. The resulting divisions
    are the closest to cubes.

    Parameters
    ----------
    dim: array_like of length 3
        Dimensions of the object to be divided.
    target_elems: int,
        Total number of elements as target for the procedure. Actual final number is
        likely to differ.
    min_val: int
        Minimum value of the number of divisions per dimension.
    strict_max: bool
        If `True`, the `target_elem` value becomes a strict maximum and the product of
        the resulting triple will be strictly smaller than the target.
    parity: {None, 'odd', 'even'}
        All elements of the resulting triple will match the given parity. If `None`, no
        parity check is performed.

    Returns
    -------
    numpy.ndarray of length 3
        array corresponding of the number of divisions for each dimension

    Examples
    --------
    >>> cells_from_dimension([1, 2, 6], 926, parity=None, strict_max=True)
    [ 4  9 25]  # Actual total: 900
    >>> cells_from_dimension([1, 2, 6], 926, parity=None, strict_max=False)
    [ 4  9 26]  # Actual total: 936
    >>> cells_from_dimension([1, 2, 6], 926, parity='odd', strict_max=True)
    [ 3 11 27]  # Actual total: 891
    >>> cells_from_dimension([1, 2, 6], 926, parity='odd', strict_max=False)
    [ 5  7 27]  # Actual total: 945
    >>> cells_from_dimension([1, 2, 6], 926, parity='even', strict_max=True)
    [ 4  8 26]  # Actual total: 832
    >>> cells_from_dimension([1, 2, 6], 926, parity='even', strict_max=False)
    [ 4 10 24]  # Actual total: 960
    """
    elems = np.prod(target_elems)  # in case target_elems is an iterable

    # define parity functions
    if parity == "odd":
        funcs = [
            lambda x, add=add, fn=fn: int(2 * fn(x / 2) + add)
            for add in (-1, 1)
            for fn in (np.ceil, np.floor)
        ]
    elif parity == "even":
        funcs = [lambda x, fn=fn: int(2 * fn(x / 2)) for fn in (np.ceil, np.floor)]
    else:
        funcs = [np.ceil, np.floor]

    # make sure the number of elements is sufficient
    elems = max(min_val**3, elems)

    # float estimate of the elements while product=target_elems and proportions are kept
    x, y, z = np.abs(dim)
    a = x ** (2 / 3) * (elems / y / z) ** (1 / 3)
    b = y ** (2 / 3) * (elems / x / z) ** (1 / 3)
    c = z ** (2 / 3) * (elems / x / y) ** (1 / 3)
    a, b, c = apportion_triple((a, b, c), min_val=min_val)
    epsilon = elems
    # run all combinations of rounding methods, including parity matching to find the
    # closest triple with the target_elems constrain
    result = [funcs[0](k) for k in (a, b, c)]  # first guess
    for funcs in product(*[funcs] * 3):
        res = [f(k) for f, k in zip(funcs, (a, b, c))]
        epsilon_new = elems - np.prod(res)
        if np.abs(epsilon_new) <= epsilon and all(r >= min_val for r in res):
            if not strict_max or epsilon_new >= 0:
                epsilon = np.abs(epsilon_new)
                result = res
    return np.array(result).astype(int)

def mesh_cylinder(object):

    n = object.meshing
    
    if isinstance(object, CylinderSegment):
        r1, r2, h, phi1, phi2 = object.dimension
    elif isinstance(object, Cylinder):
        r1, r2, h, phi1, phi2 = (
            0,
            object.dimension[0] / 2,
            object.dimension[1],
            0,
            360,
        )
    else:
        raise TypeError("Input must be a Cylinder or CylinderSegment")

    al = (r2 + r1) * 3.14 * (phi2 - phi1) / 360  # arclen = D*pi*arcratio
    dim = al, r2 - r1, h
    # "unroll" the cylinder and distribute the target number of elemens along the
    # circumference, radius and height.
    nphi, nr, nh = cells_from_dimension(dim, n)

    elems = np.prod([nphi, nr, nh])

    r = np.linspace(r1, r2, nr + 1)
    dh = h / nh
    cells = []
    for r_ind in range(nr):
        # redistribute number divisions proportionally to the radius
        nphi_r = max(1, int(r[r_ind + 1] / ((r1 + r2) / 2) * nphi))
        phi = np.linspace(phi1, phi2, nphi_r + 1)
        for h_ind in range(nh):
            pos_h = dh * h_ind - h / 2 + dh / 2
            # use a cylinder for the innermost cells if there are at least 3 layers and
            # if it is closed, use cylinder segments otherwise
            if nr >= 3 and r[r_ind] == 0 and phi2 - phi1 == 360:
                cell = (0,0,pos_h)
                cells.append(cell)
            else:
                for phi_ind in range(nphi_r):
                    dimension = (
                        r[r_ind],
                        r[r_ind + 1],
                        dh,
                        phi[phi_ind],
                        phi[phi_ind + 1],
                    )
                    radial_coord = (r[r_ind] + r[r_ind + 1]) / 2
                    angle_coord = (phi[phi_ind] + phi[phi_ind + 1]) / 2
                    
                    cell = (radial_coord * np.cos(np.deg2rad(angle_coord)), 
                            radial_coord * np.sin(np.deg2rad(angle_coord)), 
                            pos_h)
                    cells.append(cell)
    #return _collection_from_obj_and_cells(cylinder, cells, **kwargs)
    return np.array(cells)


# def mesh_cylinder(object):
#     """
#     create cylinder mesh from object meshing parameter
#     """
#     n = object.meshing
#     dia, height = object.dimension

#     a1 = -dia/2+dia/(2*n)
#     b1 =  dia/2-dia/(2*n)
#     a2 = -height/2+height/(2*n)
#     b2 =  height/2-height/(2*n)

#     if dia > height:
#         c1 = n
#         c2 = int(n/dia*height)
#         if c2<3:
#             print("Warning: bad cylinder mesh ratio. Increase meshing parameter.")
#             c2 = 3
#     else:
#         c2 = n
#         c1 = int(n/height*dia)
#         if c1<3:
#             print("Warning: bad cylinder mesh ratio. Increase meshing parameter.")
#             c1 = 3

#     mesh = np.mgrid[a1:b1:c1*1j, a1:b1:c1*1j, a2:b2:c2*1j].T.reshape(c1*c1*c2,3)
#     mask = np.linalg.norm(mesh[:,:2], axis=1) < dia/2
#     return mesh[mask]


# def mesh_cylinder_segment(object):
#     """
#     create cylinder mesh from object meshing parameter
#     """
#     n = object.meshing
#     r1, r2, h, phi1, phi2 = object.dimension

#     a1 = -r2+r2/n
#     b1 =  r2-r2/n
#     a2 = -h/2 + h/(2*n)
#     b2 =  h/2 - h/(2*n)

#     dia = 2*r2
#     if dia > h:
#         c1 = n
#         c2 = int(n/dia*h)
#         c2 = 3 if c2<3 else c2
#     else:
#         c2 = n
#         c1 = int(n/h*dia)
#         c1 = 3 if c1<3 else c1

#     mesh = np.mgrid[a1:b1:c1*1j, a1:b1:c1*1j, a2:b2:c2*1j].T.reshape(c1*c1*c2,3)
#     import magpylib as magpy

#     mask1 = np.linalg.norm(mesh[:,:2], axis=1) > r1
#     mask2 = np.linalg.norm(mesh[:,:2], axis=1) < r2
#     mask3 = np.arctan2(mesh[:,0], mesh[:,1]) > phi1/180*np.pi
#     mask4 = np.arctan2(mesh[:,0], mesh[:,1]) < phi2/180*np.pi
#     return mesh[mask1*mask2*mask3*mask4]


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
