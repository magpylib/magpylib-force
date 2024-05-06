# Magpylib-Force

This package provides force and torque computation for Magpylib objects. Its work in progress and will finally be integrate into the Magpylib package.

# Installation

`pip install magpylib-force`

# Example code:

```
import numpy as np
import magpylib as magpy
from magpylib_force import getFT

source = magpy.magnet.Sphere(diameter=1, polarization=(1,2,3))

line = magpy.current.Polyline(
    current=1,
    vertices=((0,0,3), (-2,0,-1), (2,0,-1), (0,0,3)),
)
line.meshing = 10 # mesh points in each segment

cube = magpy.magnet.Cuboid(
    dimension=(1,1,1),
    polarization=(0,0,1),
    position=(1,2,3)
)
cube.meshing = (2,2,2) # mesh grid

FT = getFT(source, [line, cube])
print(FT)
```
