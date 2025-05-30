from __future__ import annotations

import warnings

import magpylib as magpy
import numpy as np

from magpylib_force.force import getFT


def test_ANSYS_cube_cube():
    """
    compare to ANSYS values using EXP04
    """
    dat = np.array(
        (
            (
                1,
                7,
                2,
                377.854600016549,
                -188.542528864873,
                -92.8990855112117,
                -313.999320817873,
                141.297867367292,
                -298.077476393508,
                2.79328118297203,
            ),
            (
                1,
                7,
                4,
                222.97324592523,
                -196.942454862682,
                -45.7293146685391,
                -94.0189749520763,
                103.019653363118,
                -499.3139545593,
                14.0657284644094,
            ),
            (
                1,
                10,
                2,
                109.392835475382,
                -40.734672622913,
                -19.9251847187105,
                -99.5513230126722,
                43.6991462807681,
                -97.0400812701753,
                0.963515112458678,
            ),
            (
                1,
                10,
                4,
                82.5409955281411,
                -57.0911591058395,
                -13.7888753933175,
                -57.9955378465955,
                37.2506186840298,
                -167.207747894386,
                1.88952437302015,
            ),
            (
                3,
                7,
                2,
                268.044882092576,
                -122.351045384398,
                -181.533246425662,
                -154.673724179409,
                370.592326952291,
                -258.784762668333,
                4.13894019162516,
            ),
            (
                3,
                7,
                4,
                164.713887772005,
                -128.609022516067,
                -95.5090438696884,
                -38.319795035053,
                268.387524559781,
                -363.857009255895,
                3.85642340519037,
            ),
            (
                3,
                10,
                2,
                87.2315925006503,
                -30.1457939660204,
                -48.080701791827,
                -66.2482297994683,
                142.329658703303,
                -76.1110023199789,
                -5.72426828151665,
            ),
            (
                3,
                10,
                4,
                66.3151258655567,
                -44.4688881017354,
                -32.4675166815565,
                -36.9604419628989,
                104.515519563477,
                -150.862082328095,
                3.51331256342416,
            ),
        )
    )
    dat[:, (0, 1, 2)] = dat[:, (2, 0, 1)]  # correct bad xyz-order

    tgt_pos = dat[:, (0, 1, 2)]
    F_fe = dat[:, (4, 5, 6)]
    T_fe = dat[:, (7, 8, 9)]
    gen = magpy.magnet.Cuboid(
        polarization=(0, 0, 1),
        dimension=(5, 5, 5),
    )
    tgt = magpy.magnet.Cuboid(
        dimension=(0.3, 0.3, 0.3),
        polarization=(0, 0, 1),
    )
    tgt.meshing = (10, 10, 10)

    for i, poz in enumerate(tgt_pos):
        tgt.position = poz
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            F, T = getFT(gen, tgt)
        T *= -1  # bad sign in original implementation

        errF = np.linalg.norm(F - F_fe[i]) / np.linalg.norm(F)
        assert errF < 0.04
        errT = np.linalg.norm(T - T_fe[i]) / np.linalg.norm(T)
        assert errT < 0.25


def test_ANSYS_loop_loop():
    """
    compare to ANSYS loop-loop computation
    Warning: ANSYS computes force very inaccurately and torque is completely off
    """
    data = (
        (
            0.5,
            0.25,
            -500,
            -0.378990675610971,
            1.16471535508468,
            -3.93840492731175,
            0.316391999895969,
            -0.323117059194582,
            0.554105980303591,
        ),
        (
            0.5,
            0.25,
            500,
            0.52319608703649,
            0.477119683892785,
            -3.95525312585684,
            -0.344823150943224,
            -0.32983037884051,
            0.548719111194436,
        ),
        (
            0.5,
            1.75,
            -500,
            -0.247318906984941,
            -0.475750180468637,
            2.20298358334411,
            -0.112567180482201,
            0.0876390943086132,
            0.40936448187346,
        ),
        (
            0.5,
            1.75,
            500,
            0.271711794735348,
            -0.229983825017913,
            2.51566260990816,
            0.091544957851064,
            0.0899073023482688,
            0.413415626972563,
        ),
        (
            1.5,
            0.25,
            -500,
            0.938566291200263,
            -2.02436288094611,
            -7.4080821527953,
            -0.618507366424736,
            1.99422271236167,
            3.65138404864739,
        ),
        (
            1.5,
            0.25,
            500,
            -0.00216325664572926,
            -2.1084068976844,
            -6.54170362000025,
            0.5995369738223,
            2.00706353179035,
            3.67068375878204,
        ),
        (
            1.5,
            0.5,
            500,
            -0.523860907171043,
            -1.99575639650359,
            -0.940337877715457,
            0.429197197716147,
            1.28142571497163,
            1.8011336159446,
        ),
        (
            1.5,
            1.75,
            -500,
            -0.328463925710092,
            -0.738926572407658,
            3.6506362987988,
            -0.0893135956299432,
            0.26486293279763,
            0.263853801870998,
        ),
        (
            1.5,
            1.75,
            500,
            -0.143145542837849,
            -0.277054173217136,
            2.16655477926914,
            0.0843493800576269,
            0.246120004205157,
            0.234717784641934,
        ),
    )

    i_squ = 0.5
    i_circ = 10

    verts1 = (
        np.array(
            (
                (0.5, 0.5, 0),
                (-0.5, 0.5, 0),
                (-0.5, -0.5, 0),
                (0.5, -0.5, 0),
                (0.5, 0.5, 0),
            )
        )
        * 1e-3
    )
    sloop = magpy.current.Polyline(
        vertices=verts1,
        current=i_squ,
    )
    sloop.meshing = 100
    ts = np.linspace(0, 2 * np.pi, 100)
    verts2 = 1.975 * np.array([(np.cos(t), np.sin(t), 0) for t in ts]) * 1e-3
    cloop = magpy.current.Polyline(vertices=verts2, current=i_circ)
    cloop.meshing = 3

    for d in data:
        c1y, c1z, c1x = d[:3]
        pos = np.array((c1x * 1e-3, c1y, c1z)) * 1e-3
        cloop.position = pos

        # fem force
        F2 = d[6:9]

        # analytical force
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            F3, _ = getFT(sources=cloop, targets=sloop)
        F3 *= 1e6

        err = np.linalg.norm(F2 - F3) / np.linalg.norm(F3)
        assert err < 0.2


def test_ANSYS_loop_magnet():
    """
    compare to FEM solution
    """
    # "yy [mm]","zz [mm]","xx [um]","F_magnet.Force_mag [mNewton]","F_magnet.Force_x [mNewton]","F_magnet.Force_y [mNewton]","F_magnet.Force_z [mNewton]","F_square.Force_mag [mNewton]","F_square.Force_x [mNewton]","F_square.Force_y [mNewton]","F_square.Force_z [mNewton]","Fv_magnet.Force_mag [mNewton]","Fv_magnet.Force_x [mNewton]","Fv_magnet.Force_y [mNewton]","Fv_magnet.Force_z [mNewton]"
    dataF = np.array(
        (
            (
                0.8,
                0.1,
                200,
                6.73621798696924e-15,
                -3.584546668289e-15,
                4.7801359239701e-16,
                -5.68323507839592e-15,
                16.6586082371456,
                -11.831409840104,
                4.03959186599999,
                11.0094807847751,
                16.6805051117026,
                11.8433422379711,
                -4.05649710352454,
                -11.0235804829885,
            ),
            (
                0.8,
                0.1,
                800,
                1.92251067162549e-15,
                4.23204435527338e-16,
                1.17475208106657e-15,
                -1.46181491177701e-15,
                13.3809780256917,
                4.44706302385542,
                5.81654806197601,
                11.2000880366462,
                13.4481224103919,
                -4.41614382881405,
                -5.86900194663829,
                -11.2651891328316,
            ),
            (
                1,
                0.1,
                800,
                0,
                0,
                0,
                0,
                10.5206977015748,
                3.18905970265248,
                6.23300504192764,
                7.85268275738579,
                10.5186835277485,
                -3.19801296156644,
                -6.2332692626082,
                -7.84613092896142,
            ),
            (
                2,
                0.1,
                800,
                0,
                0,
                0,
                0,
                0.85304709445094,
                -0.123455805806032,
                0.800549927287909,
                -0.267521631430614,
                0.821752082865582,
                0.0413733953267691,
                -0.783525909803599,
                0.244237336456775,
            ),
            (
                0.8,
                0.8,
                200,
                0,
                0,
                0,
                0,
                2.84620865578413,
                -2.56083873331669,
                0.556365961005322,
                1.11061497002341,
                2.86935053705004,
                2.60797676976126,
                -0.418458292630111,
                -1.12094706841316,
            ),
            (
                0.8,
                0.8,
                800,
                0,
                0,
                0,
                0,
                2.54654428882465,
                -0.0827582901825677,
                0.976845061909694,
                2.35027926114625,
                2.56309182955053,
                0.160281418778475,
                -1.03704319103324,
                -2.33843772921894,
            ),
        )
    )
    # "yy [mm]","zz [mm]","xx [um]","Tvx_m.Torque [uNewtonMeter]","Tvy_m.Torque [uNewtonMeter]","Tvz_m.Torque [uNewtonMeter]","Tx_square.Torque [uNewtonMeter]","Ty_square.Torque [uNewtonMeter]","Tz_square.Torque [uNewtonMeter]"
    dataT = np.array(
        (
            (
                0.8,
                0.1,
                200,
                -0.152752678323172,
                2.53687184879894,
                0.31545404743791,
                0.140313023228343,
                -2.53736199874947,
                -0.30517463426644,
            ),
            (
                0.8,
                0.1,
                800,
                0.959480198170481,
                2.89789848899664,
                -0.435258330080741,
                -0.870512145889583,
                -2.87366675400084,
                0.430657575370342,
            ),
            (
                1,
                0.1,
                800,
                0.800553299865047,
                2.17870450804352,
                -0.381791156481078,
                -0.81345605829445,
                -2.17610385228856,
                0.398898899795887,
            ),
            (
                2,
                0.1,
                800,
                0.853704644034881,
                0.427380408643971,
                -0.0124517565704488,
                -0.888547449837403,
                -0.41747637388676,
                0.0181451681407211,
            ),
            (
                0.8,
                0.8,
                200,
                0.0940961486719588,
                1.24970685094319,
                0.00346834797401525,
                -0.118237549519446,
                -1.18243794063886,
                -0.0161649360104749,
            ),
            (
                0.8,
                0.8,
                800,
                0.340487764979492,
                0.756762736556279,
                -0.0413327662068715,
                -0.346866763509663,
                -0.707475285454173,
                0.00988830084537582,
            ),
        )
    )

    verts1 = (
        np.array(
            (
                (0.5, 0.5, 0),
                (-0.5, 0.5, 0),
                (-0.5, -0.5, 0),
                (0.5, -0.5, 0),
                (0.5, 0.5, 0),
            )
        )
        * 1e-3
    )
    loop = magpy.current.Polyline(
        vertices=verts1,
        current=50,
    )
    magnet = magpy.magnet.Cuboid(
        dimension=np.array((1, 2, 1)) * 1e-3,
        polarization=(1, 0, 0),
    )

    for dat, dat2 in zip(dataF, dataT, strict=False):
        c1y, c1z, c1x = dat[:3]
        pos = np.array((c1x * 1e-3, c1y, c1z)) * 1e-3
        pos[2] += 0.5 * 1e-3
        magnet.position = pos

        c1yb, c1zb, c1xb = dat2[:3]
        pos2 = np.array((c1xb * 1e-3, c1yb, c1zb)) * 1e-3
        pos2[2] += 0.5 * 1e-3
        np.testing.assert_allclose(pos, pos2)

        F2 = dat[8:11]
        F1 = dat[12:15]
        T2 = dat2[3:6]
        T1 = dat2[6:9]

        loop.meshing = 1000
        magnet.meshing = (10, 20, 10)
        F3, T3 = getFT(sources=loop, targets=magnet, anchor=(0, 0, 0))
        T3 *= -1  # bad sign at initial test design
        F4, T4 = getFT(sources=magnet, targets=loop, anchor=(0, 0, 0))
        T4 *= -1  # bad sign at initial test design
        F3 *= 1e3
        F4 *= 1e3
        T3 *= 1e6
        T4 *= 1e6

        err = np.linalg.norm(F2 + F3) / np.linalg.norm(F3)
        assert err < 0.15
        err = np.linalg.norm(F1 + F4) / np.linalg.norm(F4)
        assert err < 0.15
        err = np.linalg.norm(T2 + T3) / np.linalg.norm(T3)
        assert err < 0.15
        err = np.linalg.norm(T1 + T4) / np.linalg.norm(T4)
        assert err < 0.15


def test_ANSYS_magnet_current_close():
    """current loop close to magnet"""

    magnet = magpy.magnet.Cuboid(
        dimension=np.array((0.5, 10, 0.3)) * 1e-3,
        polarization=(1, 0, 0),
    )
    magnet.meshing = (5, 50, 3)

    # wire spit up into 4 parts
    d = 0.025  # put PolyLine in center of crosssection
    t = 0.01  # put PolyLine in center of crosssection
    verts1 = (
        np.array(
            (
                (-0.25 + d, -4 + d, t),
                (0.25 - d, -4 + d, t),
                (0.25 - d, 4 - d, t),
                (-0.25 + d, 4 - d, t),
                (-0.25 + d, -4 + d, t),
            )
        )
        * 1e-3
    )
    discr = 10 * 1e3  # wire discretizations per meter
    wires = []
    for i in range(4):
        wire = magpy.current.Polyline(vertices=(verts1[i : i + 2]))
        mw = int(discr * np.linalg.norm(verts1[i] - verts1[i + 1])) + 1
        wire.meshing = mw
        wires.append(wire)

    # "I_square [mA]","yy [mm]","zz [mm]","xx [um]","F_square.Force_x [uNewton]","F_square.Force_y [uNewton]","F_square.Force_z [uNewton]","Fv_magnet.Force_x [uNewton]","Fv_magnet.Force_y [uNewton]","Fv_magnet.Force_z [uNewton]"
    datF = np.array(
        (
            (
                50,
                0,
                0.2,
                200,
                -38.697482833418,
                0.0002242350531848,
                59.4310879419995,
                29.0564031184439,
                7.11803882106958,
                -47.9368393911027,
            ),
            (
                50,
                0,
                0.2,
                500,
                31.1552670207915,
                -0.00142810245237877,
                29.8624423913996,
                -65.3924604890968,
                11.4789019916509,
                -51.8408330694864,
            ),
            (
                50,
                0,
                0.2,
                800,
                15.1908535165256,
                -3.87787107323596e-05,
                -2.59112312011276,
                -10.2503856180475,
                3.77071833140656,
                -40.2905673214654,
            ),
            (
                50,
                1,
                0.2,
                200,
                -38.2166912173904,
                0.74775728098968,
                59.1285115576116,
                -7.1275613537107,
                -8.16434162283945,
                -70.8158765538541,
            ),
            (
                50,
                1,
                0.2,
                500,
                30.9248391472397,
                0.711934423270261,
                29.7711197385536,
                -47.0335059106064,
                5.30602983129916,
                -94.6243138266417,
            ),
            (
                50,
                1,
                0.2,
                800,
                14.9818522746782,
                0.285580482211302,
                -2.48002525728516,
                -17.4778153030685,
                14.3839498468311,
                -33.3455008748144,
            ),
            (
                50,
                0,
                0.5,
                200,
                -15.4440842393698,
                0.000630278271938781,
                13.8040620063476,
                -29.3586587495271,
                -7.04552102088018,
                -31.9675710200005,
            ),
            (
                50,
                0,
                0.5,
                500,
                2.40460818037774,
                -0.000453637988892101,
                14.3771835667471,
                59.7897254642478,
                -7.72738474941987,
                -5.30574925396906,
            ),
            (
                50,
                0,
                0.5,
                800,
                6.6367471620929,
                0.000255164078938144,
                4.82949919020895,
                -18.5474946200134,
                1.33369577647253,
                -10.881515270821,
            ),
            (
                50,
                1,
                0.5,
                200,
                -15.145226144758,
                0.235812935580012,
                13.6277157811128,
                -1.85120493391665,
                -2.66991565009376,
                -10.0791303759911,
            ),
            (
                50,
                1,
                0.5,
                500,
                2.40322610295551,
                0.343272487193449,
                14.182570171465,
                42.5947990318086,
                8.34168525083406,
                -3.27399171150132,
            ),
            (
                50,
                1,
                0.5,
                800,
                6.52589219503979,
                0.243370384679762,
                4.7572297397883,
                -42.1269990492652,
                -7.52176458860419,
                -4.45812096051993,
            ),
        )
    )
    # "I_square [mA]","yy [mm]","zz [mm]","xx [um]","Tvx_m.Torque [nNewtonMeter]","Tvy_m.Torque [nNewtonMeter]","Tvz_m.Torque [nNewtonMeter]","Tx_square.Torque [nNewtonMeter]","Ty_square.Torque [nNewtonMeter]","Tz_square.Torque [nNewtonMeter]"
    datT = np.array(
        (
            (
                50,
                0,
                0.2,
                200,
                100.376931034983,
                1.29132310544713,
                -44.5742433626561,
                0.0135554469175028,
                -12.259497183979,
                0.0234822598655818,
            ),
            (
                50,
                0,
                0.2,
                500,
                -35.3674218480089,
                15.2488447927969,
                16.5356714838072,
                0.0128574072989242,
                -2.11323611976223,
                0.0140856152834011,
            ),
            (
                50,
                0,
                0.2,
                800,
                78.1145191288726,
                28.9160470934922,
                42.7406950893594,
                -0.00270261511070544,
                3.82101753288969,
                0.001268285244331,
            ),
            (
                50,
                1,
                0.2,
                200,
                253.962525714332,
                -3.58528490354542,
                -18.3631031181225,
                0.906354335079891,
                -11.983108298508,
                1.88954778000799,
            ),
            (
                50,
                1,
                0.2,
                500,
                4.46861420652312,
                17.1480742529462,
                29.6735385079683,
                0.0283925165409689,
                -2.05547893566587,
                -0.837502175937236,
            ),
            (
                50,
                1,
                0.2,
                800,
                -134.966632069798,
                26.1191706299256,
                -3.58074585018869,
                -0.728913411271276,
                3.76373308619741,
                -0.751039454165596,
            ),
            (
                50,
                0,
                0.5,
                200,
                -23.9024071346421,
                -17.6356253739313,
                -183.322034423475,
                0.00115987401080014,
                -6.23683839052518,
                0.00221598446713591,
            ),
            (
                50,
                0,
                0.5,
                500,
                2.78526886521961,
                31.7179118567721,
                -24.2255502578008,
                0.00550707289213479,
                -2.67769582918341,
                -0.00468607713562535,
            ),
            (
                50,
                0,
                0.5,
                800,
                -40.0407569434816,
                3.37798020460892,
                -39.931893857586,
                -0.00396662587435437,
                0.350027477836724,
                -0.00218198771623613,
            ),
            (
                50,
                1,
                0.5,
                200,
                -2.0307023796287,
                0.583836236430798,
                75.3662429442374,
                0.542371195596806,
                -6.0801890938976,
                1.19222997481852,
            ),
            (
                50,
                1,
                0.5,
                500,
                23.698487998103,
                18.7804224066191,
                51.5633199828762,
                0.51802382216635,
                -2.60537736999987,
                0.0549590129589455,
            ),
            (
                50,
                1,
                0.5,
                800,
                145.274705178714,
                -19.4722759246531,
                6.63502773049354,
                0.0223391263649002,
                0.348069100595905,
                -0.392873453263144,
            ),
        )
    )

    for d, t in zip(datF, datT, strict=False):
        i0 = d[0] * 1e-3  # ampere
        pos = np.array((d[3] * 1e-3, d[1], d[2])) * 1e-3
        f2 = np.array((d[4], d[5], d[6])) * 1e-6
        # f1 = np.array((d[7], d[8], d[9])) * 1e-6 # TODO check if necessary

        # t1 = np.array((t[4], t[5], t[6])) * 1e-9 # TODO check if necessary
        t2 = np.array((t[7], t[8], t[9])) * 1e-9

        for wire in wires:
            wire.current = i0
        magnet.position = pos + np.array((0, 0, 0.15)) * 1e-3

        F1, _ = getFT(wires, magnet, anchor=(0, 0, 0))
        F2, T2 = np.sum(getFT(magnet, wires, anchor=(0, 0, 0)), axis=0)
        T2 *= -1  # bad sign at initial test design

        assert np.linalg.norm(F1 + F2) / np.linalg.norm(F1) < 1e-3
        assert np.linalg.norm(f2 - F2) / np.linalg.norm(F2) < 1e-2
        assert np.linalg.norm(t2 + T2) / np.linalg.norm(T2) < 0.1
