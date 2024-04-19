import numpy as np
import magpylib as magpy
from mforce.force import getFTcube


def test_ANSYS_cube_cube():
    """
    compare to ANSYS values using EXP04
    """
    dat = np.array((
        (1,7,2,377.854600016549,-188.542528864873,-92.8990855112117,-313.999320817873,141.297867367292,-298.077476393508,2.79328118297203),
        (1,7,4,222.97324592523,-196.942454862682,-45.7293146685391,-94.0189749520763,103.019653363118,-499.3139545593,14.0657284644094),
        (1,10,2,109.392835475382,-40.734672622913,-19.9251847187105,-99.5513230126722,43.6991462807681,-97.0400812701753,0.963515112458678),
        (1,10,4,82.5409955281411,-57.0911591058395,-13.7888753933175,-57.9955378465955,37.2506186840298,-167.207747894386,1.88952437302015),
        (3,7,2,268.044882092576,-122.351045384398,-181.533246425662,-154.673724179409,370.592326952291,-258.784762668333,4.13894019162516),
        (3,7,4,164.713887772005,-128.609022516067,-95.5090438696884,-38.319795035053,268.387524559781,-363.857009255895,3.85642340519037),
        (3,10,2,87.2315925006503,-30.1457939660204,-48.080701791827,-66.2482297994683,142.329658703303,-76.1110023199789,-5.72426828151665),
        (3,10,4,66.3151258655567,-44.4688881017354,-32.4675166815565,-36.9604419628989,104.515519563477,-150.862082328095,3.51331256342416),
    ))
    dat[:, (0,1,2)] = dat[:, (2,0,1)] # correct bad xyz-order

    tgt_pos = dat[:,(0,1,2)]
    F_fe = dat[:,(4,5,6)]
    T_fe = dat[:,(7,8,9)]
    gen = magpy.magnet.Cuboid(
        polarization=(0,0,1),
        dimension=(5,5,5),
    )
    tgt = magpy.magnet.Cuboid(
        dimension=(.3,.3,.3),
        polarization=(0,0,1),
    )
    tgt.mesh=(10,10,10)

    for i,poz in enumerate(tgt_pos):
        tgt.position = poz
        F,T = getFTcube(gen, tgt)

        errF = np.linalg.norm(F - F_fe[i])/np.linalg.norm(F)
        assert errF < 0.04
        errT = np.linalg.norm(T - T_fe[i])/np.linalg.norm(T)
        assert errT < 0.25