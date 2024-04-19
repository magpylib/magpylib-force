# Magpylib-Force

This package provides force and torque computation for Magpylib objects. Its work in progress and will finally be integrate into the Magpylib package.

# Installation

`pip install magpylib-force`

# Example code:

```
import magpylib as magpy
from mforce import getFTcube

# all inputs and output in SI units

cube1 = magpy.magnet.Cuboid(
    dimension=(1,1,1),
    polarization=(1,1,1),
)
cube2 = cube1.copy(position=(1,0,2))
cube2.mesh = (10,10,10)

F,T = getFTcube(cube1, [cube2, cube2.copy()], anchor=(0,0,0))
print(F)
print(T)
```
