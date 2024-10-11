"""
Force computation codes
"""

import magpylib as magpy
import numpy as np

from magpylib._src.obj_classes.class_current_Polyline import Polyline
from magpylib._src.obj_classes.class_magnet_Cuboid import Cuboid
from magpylib._src.obj_classes.class_magnet_Sphere import Sphere
from magpylib._src.obj_classes.class_magnet_Cylinder import Cylinder
from magpylib._src.obj_classes.class_magnet_CylinderSegment import CylinderSegment

from magpylib_force.meshing import mesh_target
from magpylib_force.utility import check_input_anchor
from magpylib_force.utility import check_input_targets


#pylint: disable=invalid-name
#pylint: disable=too-many-locals


def getFT(sources, targets, anchor=None, eps=1e-5, squeeze=True):
    """
    Compute magnetic force and torque acting on the targets that are exposed
    to the magnetic field of the sources.
    SI units are assumed for all inputs and outputs.

    Parameters
    ----------
    sources: source and collection objects or 1D list thereof
        Sources that generate the magnetic field. Can be a single source (or collection)
        or a 1D list of l sources and/or collection objects.

    targets: single object or 1D list of t objects that are Sphere, Cuboid, Polyline,
        Cylinder, or CylinderSegment. Force and Torque acting on targets in the magnetic
        field generated by the sources will be computed. A target must have a valid
        `meshing` parameter.

    anchor: array_like, default=None
        The Force adds to the Torque via the anchor point. For a freely floating magnet
        this would be the barycenter. If `anchor=None`, this part of the Torque computation
        is ommitted and a warning is thrown.

    eps: float, default=1e-5
        This is only used for magnet targets for computing the magnetic field gradient
        using finite differences (FD). `eps` is the FD step size. A good value
        is 1e-5 * characteristic_system_size (magnet size, distance between sources
        and targets, ...).

    squeeze: bool, default=True
        The output of the computation has the shape (n,3) where n corresponds to the number
        of targets. By default this is reduced to (3,) when there is only one target.

    Returns
    -------
    Force-Torque as ndarray of shape (2,3), or (t,2,3) when t targets are given
    """
    anchor = check_input_anchor(anchor)
    targets = check_input_targets(targets)
    # MISSING: allow Collections as targets

    n = len(targets)

    # split targets into lists of similar types
    TARGET_TYPES = [Cuboid, Polyline, Sphere, Cylinder, CylinderSegment]
    getFT_FUNCS = [getFTmagnet, getFTcurrent, getFTmagnet, getFTmagnet, getFTmagnet]
    objects = [[] for _ in TARGET_TYPES]
    orders  = [[] for _ in TARGET_TYPES]

    for i,tgt in enumerate(targets):
        for j,ttyp in enumerate(TARGET_TYPES):
            if isinstance(tgt, ttyp):
                objects[j].append(tgt)
                orders[j].append(i)

    # allocate FT
    FT = np.zeros((n, 2, 3))

    # FT-computation and broadcasting
    for i in range(len(TARGET_TYPES)):
        if objects[i]:
            ft_part = getFT_FUNCS[i](sources, objects[i], eps=eps, anchor=anchor)
            ft_part = np.swapaxes(ft_part, 0, 1)
            for ft,j in zip(ft_part, orders[i]):
                FT[j] = ft

    if squeeze:
        return np.squeeze(FT)
    return FT


def volume(target):
    """
    utility for getFTmagnet: compute physical volume of target magnets
    target: target magnet object
    return: volume, float
    """
    if isinstance(target, Cuboid):
        return np.prod(target.dimension)
    if isinstance(target, Sphere):
        return target.diameter**3 * np.pi / 6
    if isinstance(target, Cylinder):
        d, h = target.dimension
        return d**2 * np.pi * h / 4
    if isinstance(target, CylinderSegment):
        r1,r2,h,phi1,phi2 = target.dimension
        return (r2**2-r1**2)*np.pi*h * (phi2-phi1)/360
    raise RuntimeError("fktn `volume` - I shouldt be here.")


def getFTmagnet(sources, targets, eps=1e-5, anchor=None):
    """
    Compute force and torque acting on magnets of same type

    Parameters
    ----------
    sources: source and collection objects or 1D list thereof
        Sources that generate the magnetic field. Can be a single source (or collection)
        or a 1D list of l sources and/or collection objects.

    targets: Cuboid object or 1D list of Cuboid objects
        Force and Torque acting on targets in the magnetic field generated by the sources
        will be computed. A target must have a valid `meshing` parameter.

    eps: float, default=1e-5
        The magnetic field gradient is computed using finite differences (FD). eps is
        the FD step size. A good value is 1e-5 * characteristic system size (magnet size,
        distance between sources and targets, ...).

    anchor: array_like, default=None
        The Force adds to the Torque via the anchor point. For a freely floating magnet
        this would be the barycenter. If `anchor=None`, this part of the Torque computation
        is ommitted.
    """
    # number of magnets
    tgt_number = len(targets)

    # create meshes
    meshes = [mesh_target(tgt) for tgt in targets]

    # number of instances of each magnet
    inst_numbers = [len(mesh) for mesh in meshes]

    # total number of instances
    no_inst = np.sum(inst_numbers)

    # cumsum of number of instances (used for indexing)
    insti = np.r_[0, np.cumsum(inst_numbers)]

    # field computation positions (1xfor B, 6x for gradB)
    POSS = np.zeros((no_inst,7,3))

    # moment of each instance
    MOM = np.zeros((no_inst,3))

    # MISSING: eps should be defined relative to the sizes of the objects
    eps_vec = np.array([(0,0,0),(eps,0,0),(-eps,0,0),(0,eps,0),(0,-eps,0),(0,0,eps),(0,0,-eps)])

    for i,tgt in enumerate(targets):
        tgt_vol = volume(tgt)
        inst_mom = tgt.orientation.apply(tgt.magnetization) * tgt_vol / inst_numbers[i]
        MOM[insti[i]:insti[i+1]] = inst_mom

        mesh = meshes[i]
        #mesh_target(tgt)
        #import matplotlib.pyplot as plt
        #ax = plt.figure().add_subplot(projection='3d')
        #ax.plot(mesh[:,0], mesh[:,1], mesh[:,2], ls='', marker='.')

        mesh = tgt.orientation.apply(mesh)
        #ax.plot(mesh[:,0], mesh[:,1], mesh[:,2], ls='', marker='.', color='r')
        #plt.show()
        #import sys
        #sys.exit()

        for j,ev in enumerate(eps_vec):
            POSS[insti[i]:insti[i+1],j] = mesh + ev + tgt.position

    BB = magpy.getB(sources, POSS, sumup=True)
    gradB = (BB[:,1::2]-BB[:,2::2]) / (2*eps)
    gradB = np.swapaxes(gradB,0,1)

    Fs = np.sum((gradB*MOM),axis=2).T
    #Ts = np.zeros((no_inst,3))
    Ts = np.cross(BB[:,0], MOM)
    if anchor is not None:
        Ts -= np.cross(POSS[:,0]-anchor, Fs)

    T = np.array([np.sum(Ts[insti[i]:insti[i+1]],axis=0) for i in range(tgt_number)])
    F = np.array([np.sum(Fs[insti[i]:insti[i+1]],axis=0) for i in range(tgt_number)])

    return np.array((F, -T))

#pylint: disable=unused-argument
def getFTcurrent(sources, targets, anchor=None, eps=None):
    """
    compute force acting on tgt Polyline
    eps is a dummy variable that is not used

    info:
    targets = Polyline objects
    segements = linear segments within Polyline objects
    instances = computation instances, each segment is split into `meshing` points
    """

    # number of Polylines
    tgt_number = len(targets)

    # segments of each Polyline
    seg_numbers = np.array([len(tgt.vertices)-1 for tgt in targets])

    # number of mesh-points of each Polyline
    mesh_numbers = np.array([tgt.meshing for tgt in targets])

    # number of instances of each Polyline
    inst_numbers = seg_numbers*mesh_numbers

    # total number of instances
    no_inst = np.sum(inst_numbers)

    # cumsum of number of instances (used for indexing)
    insti = np.r_[0, np.cumsum(inst_numbers)]

    # path vector of each instance
    LVEC = np.zeros((no_inst,3))
    # central location of each instance
    POSS = np.zeros((no_inst,3))
    # current of each instance
    CURR = np.zeros((no_inst,))

    for i,tgt in enumerate(targets):
        verts = tgt.orientation.apply(tgt.vertices)
        mesh = mesh_numbers[i]

        lvec = np.repeat(verts[1:] - verts[:-1], mesh, axis=0)/mesh
        LVEC[insti[i]:insti[i+1]] = lvec

        CURR[insti[i]:insti[i+1]] = [tgt.current]*mesh*seg_numbers[i]

        for j in range(seg_numbers[i]):
            #pylint: disable=line-too-long
            poss = np.linspace(verts[j]+lvec[j*mesh]/2, verts[j+1]-lvec[j*mesh]/2, mesh) + tgt.position
            POSS[insti[i]+mesh*j:insti[i]+mesh*(j+1)] = poss

    # field of every instance
    B = magpy.getB(sources, POSS, sumup=True)

    # force on every instance
    F = (CURR * np.cross(LVEC, B).T).T

    # torque on every instance + sumup for every target
    if anchor is not None:
        T = np.cross(anchor - POSS, F)
        T = np.array([np.sum(T[insti[i]:insti[i+1]],axis=0) for i in range(tgt_number)])
    else:
        T = np.zeros((tgt_number,3))

    # sumup force for every target
    F = np.array([np.sum(F[insti[i]:insti[i+1]],axis=0) for i in range(tgt_number)])

    return np.array((F, -T))
