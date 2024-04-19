import numpy as np
import magpylib as magpy
from magpylib_force.force import getFTcube
from magpylib_force.force import getFTwire

def test_getFTwire_target_inputs_01():
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
    wire1.mesh=20
    wire2.mesh=20
    wire3.mesh=20

    F1,T1 = getFTwire(src1, wire1, anchor=(0,0,0))
    F2,T2 = getFTwire(src1, [wire1], anchor=(0,0,0))
    F3,T3 = getFTwire(src1, [wire1,wire2,wire3], anchor=(0,0,0))
    np.testing.assert_allclose(F1,F2)
    np.testing.assert_allclose(F1,F3[0])
    np.testing.assert_allclose(F1,F3[1])
    np.testing.assert_allclose(F1,F3[2])
    
    np.testing.assert_allclose(T1,T2)
    np.testing.assert_allclose(T1,T3[0])
    np.testing.assert_allclose(T1,T3[1])
    np.testing.assert_allclose(T1,T3[2])

    F4,T4 = getFTwire(src1, wire1, anchor=(0,0,0))
    F5,T5 = getFTwire([src1, src2], wire1, anchor=(0,0,0))
    
    np.testing.assert_allclose(F4*2,F5)
    np.testing.assert_allclose(T4*2,T5)