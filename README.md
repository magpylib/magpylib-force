# Magpylib-Force

This package provides force and torque computation for Magpylib objects. It is work in progress and will finally be integrated into the Magpylib package.

This work is supported by a NUMFOCUS grant: SDG 2024 round 1

**Warning:**
- relies on magpylib version>=5.0.1
- does not work with paths
- limited target types available

# Installation

`pip install magpylib-force`

# Example code:

```
import numpy as np
import magpylib as magpy
from magpylib_force import getFT

# Note: All inputs and outputs are in SI units.

# SOURCES: Magpylib source objects that generates the field
s1 = magpy.magnet.Sphere(diameter=1, polarization=(1,2,3))
s2 = magpy.magnet.Cuboid(dimension=(1,1,1), polarization=(1,2,3))

# TARGETS: targets that experience a force in the field of all sources
cube = magpy.magnet.Cuboid(
    dimension=(1,1,1),
    polarization=(0,0,1),
    position=(1,2,3),
)
cube.meshing = (5,5,5) # regular rectangular mesh in Cuboid

line = magpy.current.Line(
    vertices=((2,2,2), (3,3,3), (3,2,-3), (2,2,2)),
    current=1,
)
line.meshing = 15 # number of cells in each segment

# Note: use only closed current loops or be prepared to violate physics

sphere = magpy.magnet.Sphere(
    diameter=(6/np.pi)**(1/3),
    polarization=(0,0,1),
    position=(2,2,2),
)
sphere.meshing = 5 # grid finesse - the higher the finer

# COMPUTE FORCE AND TORQUE
FT = getFT([s1,s2], [cube, line, sphere], anchor=(0,0,0))

# Note: output shape is (3,2,3) from (targets, FT, xyz)
print(FT)

# Note: The magnet-numbers are large because the magnets are large :)
```
