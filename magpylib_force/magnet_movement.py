import numpy as np
import magpylib as magpy
from magpylib_force import getFT
from scipy.spatial.transform import Rotation as R

def inertia_tensor_cuboid_solid(mass, dimensions):
    dimensions_sq = dimensions**2
    return mass/12 * np.array([[dimensions_sq[1]+dimensions_sq[2],0.,0.], [0.,dimensions_sq[0]+dimensions_sq[2],0.], [0.,0.,dimensions_sq[0]+dimensions_sq[1]]])

def inertia_tensor_sphere_solid(mass, diameter):
    return mass * diameter**2 / 10 * np.identity(3)

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
    dimension2 = np.array([2,1,1])
    t1 = magpy.magnet.Cuboid(position=(-2,0,0), dimension=dimension1, polarization=(1,0,0), orientation=R.from_euler('y', -40, degrees=True))
    t1.meshing = (5,5,5)
    t2 = magpy.magnet.Cuboid(position=(2,0,0), dimension=dimension1, polarization=(1,0,0), orientation=R.from_euler('y', 40, degrees=True))
    t2.meshing = (5,5,5)


    m1 = 1
    m2 = 1
    I1 = inertia_tensor_cuboid_solid(m1, dimension1)
    I2 = inertia_tensor_cuboid_solid(m2, dimension2)

    moving_system = Moving_system([t1, t2], [], np.array([m1, m2]),  np.array([I1, I2]),  np.array([[0.,0.,0.], [0.,0.,0.]]),  np.array([[0.,0.,0.], [0.,0.,0.]]))

    #moving_system.move(0.001)


    for i in range(50):
        moving_system.move(0.001)
        p = moving_system.display()
        #p.show()
        p.off_screen = True
        p.screenshot('tmp/{:04d}.png'.format(i))


    from make_gif import make_gif
    make_gif("test", duration=50)
