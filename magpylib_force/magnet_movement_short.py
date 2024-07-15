import numpy as np
import magpylib as magpy
from magpylib_force import getFT
from scipy.spatial.transform import Rotation as R

def inverse_inertia_tensor_cuboid_solid(mass, dimensions):
    dimensions_sq = dimensions**2
    inv_tensor = 12/mass * np.array([[1/(dimensions_sq[1]+dimensions_sq[2]),0.,0.], [0.,1/(dimensions_sq[0]+dimensions_sq[2]),0.], [0.,0.,1/(dimensions_sq[0]+dimensions_sq[1])]])
    return inv_tensor

def inverse_inertia_tensor_sphere_solid(mass, diameter):
    return 10 / mass / diameter**2 * np.identity(3)

def apply_movement(targets, dt):
    """defines magnet system that is capable for moving according to force and torque
    Parameters
    ----------
    targets: magpylib collection
        Target magnets where movement is performed on
    dt: float
        finite time step for movement simulation
    """
    n_targets = len(targets)

    # calculate force and torque
    FTs = np.zeros((n_targets, 2, 3))
    for i in range(n_targets):
        # sources are all magnets instead of target
        FTs[i,:,:] = getFT(targets[:i] + targets[i+1:], [targets[i]], anchor=None)

    # simulate movement
    for i in range(n_targets):
        # calculate movement and rotation
        targets[i].velocity = targets[i].velocity + dt/targets[i].mass * FTs[i,0,:]
        targets[i].angular_velocity = targets[i].angular_velocity + dt*targets[i].orientation.apply(np.dot(targets[i].inverse_inertia_tensor, targets[i].orientation.inv().apply(FTs[i,1,:])))
        targets[i].position = targets[i].position + dt * targets[i].velocity
        targets[i].orientation = R.from_rotvec(dt*targets[i].angular_velocity)*targets[i].orientation

        print('magnet', i)
        print('position after', targets[i].position)

def display(targets):

    n_targets = len(targets)
    

    p = magpy.show(targets, backend='pyvista', return_fig=True)

    for i in range(n_targets):
        # sources are all magnets instead of target
        FTs = getFT(targets[:i] + targets[i+1:], [targets[i]], anchor=None)

        force_torque_mag = np.linalg.norm(FTs, axis=-1)
        velocities_mag = np.linalg.norm(targets[i].velocity)
        angular_velocity_mag = np.linalg.norm(targets[i].angular_velocity)

        p.add_arrows(cent=targets[i].position, direction=FTs[0,:], mag=1/force_torque_mag[0], color='g')
        p.add_arrows(cent=targets[i].position, direction=targets[i].velocity, mag=1/velocities_mag, color='b')
        p.add_arrows(cent=targets[i].position, direction=FTs[1,:], mag=1/force_torque_mag[1], color='r')
        p.add_arrows(cent=targets[i].position, direction=targets[i].angular_velocity, mag=1/angular_velocity_mag, color='m')

        p.camera.position = (0, -10, 0)

    return p

    
if __name__ == "__main__":

    # TARGETS: Magpylib target objects that move according to field
    t1 = magpy.magnet.Cuboid(position=(-2,0,0), dimension=np.array([2,1,1]), polarization=(1,0,0), orientation=R.from_euler('y', -40, degrees=True))
    t1.meshing = (20,20,20)
    t1.mass = 1
    t1.inverse_inertia_tensor = inverse_inertia_tensor_cuboid_solid(t1.mass, t1.dimension)
    t1.velocity = np.array([0.,0.,0.])
    t1.angular_velocity = np.array([0.,0.,0.])
    t2 = magpy.magnet.Sphere(position=(2,0,0), diameter=1, polarization=(1,0,0), orientation=R.from_euler('y', 40, degrees=True))
    t2.meshing = 20
    t2.mass = 1
    t2.inverse_inertia_tensor = inverse_inertia_tensor_sphere_solid(t2.mass, t2.diameter)
    t2.velocity = np.array([0.,0.,0.])
    t2.angular_velocity = np.array([0.,0.,0.])


    targets = magpy.Collection(t1, t2)
    dt = 0.001

    for i in range(200):
        apply_movement(targets, dt)
        p = display(targets)
        p.show()
