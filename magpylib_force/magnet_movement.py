import numpy as np
import magpylib as magpy
from magpylib_force import getFT
from scipy.spatial.transform import Rotation as R

def inertia_tensor_cuboid_solid(mass, dimensions):
    """calculates inertia tensor of cuboid with homogeneous mass density in coordinate system parallel to the edges

    Parameters
    ----------
    mass: float or ndarray, shape (n,)
        total mass of cuboid/s

    dimensions: ndarray, shape (3,) or (n,3)
        (x,y,z) dimensions of cuboid/s 
    """

    dimensions_sq = dimensions**2
    if np.isscalar(mass):
        assert (len(dimensions) == 3) and (len(dimensions.shape) == 1), 'since "mass" is a scalar, "dimensions" must be an ndarray of length 3'
        return mass/12 * np.array([[dimensions_sq[1]+dimensions_sq[2],0.,0.], [0.,dimensions_sq[0]+dimensions_sq[2],0.], [0.,0.,dimensions_sq[0]+dimensions_sq[1]]])
    else:
        assert (len(mass) == len(dimensions)) and (dimensions.shape[1] == 3), 'length of "mass" and "dimensions" must be the same and "dimensions" must have the shape (n,3)'
        n = len(mass)
        zero_entries = np.zeros(n*3)
        # calculate the results as 1D vector, ready to be reshaped
        result = np.concatenate(((dimensions_sq[:,1]+dimensions_sq[:,2]) * mass/12,
                                 zero_entries,
                                 (dimensions_sq[:,0]+dimensions_sq[:,2]) * mass/12,
                                 zero_entries,
                                 (dimensions_sq[:,0]+dimensions_sq[:,1]) * mass/12,
                                 )) 
        return result.reshape((n,3,3), order='F')

def inertia_tensor_sphere_solid(mass, diameter):
    """calculates inertia tensor of sphere with homogeneous mass density

    Parameters
    ----------
    mass: float or ndarray, shape (n,)
        total mass of sphere/s

    diameter: float or ndarray, shape (n,)
        diameter/s of sphere/s
    """

    if np.isscalar(mass) and np.isscalar(diameter):
        return mass * diameter**2 / 10 * np.identity(3)
    elif len(mass) == len(diameter):
        n = len(mass)
        zero_entries = np.zeros(n*3)
        matrix_entry = mass * diameter**2 / 10
        # calculate the results as 1D vector, ready to be reshaped
        result = np.concatenate((matrix_entry,
                                 zero_entries,
                                 matrix_entry,
                                 zero_entries,
                                 matrix_entry,
                                 )) 
        return result.reshape((n,3,3), order='F')

    else:
        raise ValueError('dimension mismatch, variables "mass" and "diameter" should have the same length')

class Moving_system:
    def __init__(self, targets, sources, masses, inertia_tensors, velocities, angular_velocities):
        """defines magnet system that is capable for moving according to force and torque

        Parameters
        ----------
        targets: list
            Target magnets where movement is performed on

        sources: list
            additional magnet sorces that are not moved

        masses: ndarray, shape (n,)
            masses of all target magnets

        moments_of_inertia: ndarray, shape (n,)
            moment of inertias of all target magnets - all magnets are approximated by spheres in this calculation

        inertia_tensors: ndarray, shape (n,3,3)
            inertia tensors of all target magnets

        velocities: ndarray, shape (n,3)
            initial velocities of all magnets

        angular_velocities: ndarray, shape (n,3)
            initial angular velocities of all magnets

        """

        # check is enough input is given
        self.n_targets = len(targets)
        assert self.n_targets == len(masses), "number of given targets and masses are not equal"
        assert self.n_targets == len(inertia_tensors), "number of given targets and inertia_tensors are not equal"
        assert self.n_targets == len(velocities), "number of given targets and velocities are not equal"
        assert self.n_targets == len(angular_velocities), "number of given targets and angular_velocities are not equal"
        self.targets = targets
        self.sources = sources
        self.masses = masses
        self.inverse_inertia_tensors = np.linalg.inv(inertia_tensors)
        self.velocities = velocities
        self.angular_velocities = angular_velocities

    def get_force_torque(self):
        FTs = np.zeros((self.n_targets, 2, 3))
        for i in range(len(self.targets)):
            # sources are all magnets instead of target
            FTs[i,:,:] = getFT(self.sources + self.targets[:i] + self.targets[i+1:], [self.targets[i]], anchor=None)
        return FTs

    def move(self, dt):
        FTs = self.get_force_torque()

        self.velocities = self.velocities + (dt/self.masses * FTs[:,0,:].T).T
        
        for i in range(self.n_targets):

            self.angular_velocities[i,:] = self.angular_velocities[i,:] + dt*self.targets[i].orientation.apply(np.dot(self.inverse_inertia_tensors[i,:,:], self.targets[i].orientation.inv().apply(FTs[i,1,:])))
            self.targets[i].position = self.targets[i].position + dt * self.velocities[i,:]
            self.targets[i].orientation = R.from_rotvec(dt*self.angular_velocities[i,:])*self.targets[i].orientation


    def display(self):
        FTs = self.get_force_torque()
        force_torque_mag = np.linalg.norm(FTs, axis=-1)
        velocities_mag = np.linalg.norm(self.velocities, axis=-1)
        angular_velocity_mag = np.linalg.norm(self.angular_velocities, axis=-1)
        p = magpy.show(self.sources, self.targets, backend='pyvista', return_fig=True)
        for i in range(self.n_targets):
            p.add_arrows(cent=self.targets[i].position, direction=FTs[i,0,:], mag=1/force_torque_mag[i,0], color='g')
            p.add_arrows(cent=self.targets[i].position, direction=self.velocities[i,:], mag=1/velocities_mag[i], color='b')
            p.add_arrows(cent=self.targets[i].position, direction=FTs[i,1,:], mag=1/force_torque_mag[i,1], color='r')
            p.add_arrows(cent=self.targets[i].position, direction=self.angular_velocities[i,:], mag=1/angular_velocity_mag[i], color='m')

        p.camera.position = (0, -10, 0)

        return p

if __name__ == "__main__":

    # TARGETS: Magpylib target objects that move according to field
    dimension1 = np.array([2,1,1])
    diameter2 = 1
    t1 = magpy.magnet.Cuboid(position=(-2,0,0), dimension=dimension1, polarization=(1,0,0), orientation=R.from_euler('y', -40, degrees=True))
    t1.meshing = (5,5,5)
    t2 = magpy.magnet.Sphere(position=(2,0,0), diameter=diameter2, polarization=(1,0,0), orientation=R.from_euler('y', 40, degrees=True))
    t2.meshing = 5


    m1 = 1
    m2 = 1
    I1 = inertia_tensor_cuboid_solid(m1, dimension1)
    I2 = inertia_tensor_sphere_solid(m2, diameter2)


    moving_system = Moving_system([t1, t2], [], np.array([m1, m2]),  np.array([I1, I2]),  np.array([[0.,0.,0.], [0.,0.,0.]]),  np.array([[0.,0.,0.], [0.,0.,0.]]))


    for i in range(3):
        moving_system.move(0.001)
        #p = moving_system.display()
        #p.show()



