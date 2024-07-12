import numpy as np
from magpylib_force.magnet_movement import inertia_tensor_cuboid_solid, inertia_tensor_sphere_solid

def test_vectorization_cuboid():
    """
    compare scalar to vectorized version
    """
    mass = np.array([1,2])
    dimension = np.array([[1,2,3], [4,5,6]])  

    res_vec = inertia_tensor_cuboid_solid(mass, dimension)
    res_scalar1 = inertia_tensor_cuboid_solid(mass[0], dimension[0,:])
    res_scalar2 = inertia_tensor_cuboid_solid(mass[1], dimension[1,:])

    np.testing.assert_allclose(res_vec[0,:,:], res_scalar1)
    np.testing.assert_allclose(res_vec[1,:,:], res_scalar2)


def test_vectorization_sphere():
    """
    compare scalar to vectorized version
    """
    mass = np.array([1,2])
    diameter = np.array([3,4])  

    res_vec = inertia_tensor_sphere_solid(mass, diameter)
    res_scalar1 = inertia_tensor_sphere_solid(mass[0], diameter[0])
    res_scalar2 = inertia_tensor_sphere_solid(mass[1], diameter[1])

    np.testing.assert_allclose(res_vec[0,:,:], res_scalar1)
    np.testing.assert_allclose(res_vec[1,:,:], res_scalar2)