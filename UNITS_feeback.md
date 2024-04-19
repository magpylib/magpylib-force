# Question1
cannot mix sources with units and sources without units in one script. I would have expecetd that unit-less sources are automatically given SI units.

# Problem
markers doesnt work with units :()


import numpy as np
import magpylib as magpy
import matplotlib.pyplot as plt

src = magpy.current.Polyline(
    current=(1, "A"),
    vertices=([(-1,0,0),(1,0,0)], "mm"),
)
tgt = magpy.current.Polyline(
    current=(1, "A"),
    vertices=([(-1,0,0),(1,0,0)], "mm"),
    position=((0,0,5), "mm"),
)
lvec = tgt.vertices[1]-tgt.vertices[0]

poss = tgt.position + (tgt.vertices[1]+tgt.vertices[0])/2
#B = src.getB(poss)
#F = tgt.current * np.cross(lvec, B)

magpy.show(src, tgt, markers=poss)

# Problem Custom sources

kann custom source nicht vor source mit unit definieren weil die sonst keine quantity bekommt

# Problem Matplotlib backend cannot show dimensionful cubes

