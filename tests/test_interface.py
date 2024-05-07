import numpy as np
import magpylib as magpy
from magpylib_force.force import getFT


def test_getFT_target_inputs_01():
    """
    check different target and source input formats give same
    """
    src1 = magpy.magnet.Sphere(polarization=(1,2,3), diameter=1)
    src2 = src1.copy()
    wire1 = magpy.current.Polyline(
        current=1,
        vertices=((2,0,0),(-2,2,2),(-2,-2,-2)),
    )
    wire2 = wire1.copy()
    wire3 = wire1.copy()
    wire1.meshing=20
    wire2.meshing=20
    wire3.meshing=20

    FT1 = getFT(src1, wire1, anchor=(0,0,0))
    FT2 = getFT(src1, [wire1], anchor=(0,0,0))
    FT3 = getFT(src1, [wire1,wire2,wire3], anchor=(0,0,0))

    np.testing.assert_allclose(FT1,FT2)
    np.testing.assert_allclose(FT1,FT3[0])
    np.testing.assert_allclose(FT1,FT3[1])
    np.testing.assert_allclose(FT1,FT3[2])

    FT4 = getFT(src1, wire1, anchor=(0,0,0))
    FT5 = getFT([src1, src2], wire1, anchor=(0,0,0))

    np.testing.assert_allclose(FT4*2,FT5)
    np.testing.assert_allclose(FT4*2,FT5)
