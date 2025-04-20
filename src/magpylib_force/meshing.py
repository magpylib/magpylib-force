"""
Object Meshing codes
"""

from __future__ import annotations

import itertools
from itertools import product

import numpy as np
from magpylib._src.obj_classes.class_magnet_Cuboid import Cuboid
from magpylib._src.obj_classes.class_magnet_Cylinder import Cylinder
from magpylib._src.obj_classes.class_magnet_CylinderSegment import CylinderSegment
from magpylib._src.obj_classes.class_magnet_Sphere import Sphere

# pylint: disable=too-many-locals


def mesh_target(obj):
    """
    Create a mesh for target objects based on their type.
    """
    mesh_functions = {
        Cuboid: mesh_cuboid,
        Sphere: mesh_sphere,
        Cylinder: mesh_cylinder,
        CylinderSegment: mesh_cylinder,
    }

    for obj_type, mesh_func in mesh_functions.items():
        if isinstance(obj, obj_type):
            return mesh_func(obj)

    msg = "Unsupported object type for meshing."
    raise ValueError(msg)


def mesh_sphere(obj):
    """
    create sphere mesh from object meshing parameter
    """
    n = obj.meshing  # target number of elements
    n /= np.pi / 6  # sphere volume VS cube volume ratio
    n_grid = int(n ** (1 / 3))  # splitting of cube sides

    dia = obj.diameter
    a = -dia / 2 + dia / (2 * n_grid)
    b = dia / 2 - dia / (2 * n_grid)
    c = n_grid * 1j
    mesh = np.mgrid[a:b:c, a:b:c, a:b:c]
    mesh = mesh.T.reshape(n_grid**3, 3) # pylint: disable=no-member
    return mesh[np.linalg.norm(mesh, axis=1) < dia / 2]


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
    for fn in product(*[funcs] * 3):
        res = [f(k) for f, k in zip(fn, (a, b, c))]
        epsilon_new = elems - np.prod(res)
        if (
            np.abs(epsilon_new) <= epsilon
            and all(r >= min_val for r in res)
            and (not strict_max or epsilon_new >= 0)
        ):
            epsilon = np.abs(epsilon_new)
            result = res
    return np.array(result).astype(int)


def mesh_cylinder(obj):
    """
    Mesh cylinder
    """
    n = obj.meshing

    if isinstance(obj, CylinderSegment):
        r1, r2, h, phi1, phi2 = obj.dimension
    elif isinstance(obj, Cylinder):
        r1, r2, h, phi1, phi2 = (
            0,
            obj.dimension[0] / 2,
            obj.dimension[1],
            0,
            360,
        )
    else:
        msg = "Input must be a Cylinder or CylinderSegment"
        raise TypeError(msg)

    al = (r2 + r1) * 3.14 * (phi2 - phi1) / 360  # arclen = D*pi*arcratio
    dim = al, r2 - r1, h
    # "unroll" the cylinder and distribute the target number of elements along the
    # circumference, radius and height.
    nphi, nr, nh = cells_from_dimension(dim, n)

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
                cell = (0, 0, pos_h)
                cells.append(cell)
            else:
                for phi_ind in range(nphi_r):
                    radial_coord = (r[r_ind] + r[r_ind + 1]) / 2
                    angle_coord = (phi[phi_ind] + phi[phi_ind + 1]) / 2

                    cell = (
                        radial_coord * np.cos(np.deg2rad(angle_coord)),
                        radial_coord * np.sin(np.deg2rad(angle_coord)),
                        pos_h,
                    )
                    cells.append(cell)
    # return _collection_from_obj_and_cells(cylinder, cells, **kwargs)
    return np.array(cells)


def mesh_cuboid(obj):
    """
    splits cuboid into given mesh
    returns grid positions relative to cuboid position
    """

    if np.isscalar(obj.meshing):
        n1, n2, n3 = cells_from_dimension(obj.dimension, obj.meshing)
    else:
        n1, n2, n3 = obj.meshing

    a, b, c = obj.dimension
    xs = np.linspace(-a / 2, a / 2, n1 + 1)
    ys = np.linspace(-b / 2, b / 2, n2 + 1)
    zs = np.linspace(-c / 2, c / 2, n3 + 1)

    dx = xs[1] - xs[0] if len(xs) > 1 else a
    dy = ys[1] - ys[0] if len(ys) > 1 else b
    dz = zs[1] - zs[0] if len(zs) > 1 else c

    xs_cent = xs[:-1] + dx / 2 if len(xs) > 1 else xs + dx / 2
    ys_cent = ys[:-1] + dy / 2 if len(ys) > 1 else ys + dy / 2
    zs_cent = zs[:-1] + dz / 2 if len(zs) > 1 else zs + dz / 2

    return np.array(list(itertools.product(xs_cent, ys_cent, zs_cent)))
