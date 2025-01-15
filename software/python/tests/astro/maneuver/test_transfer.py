import numpy as np
import pytest

import src.valladopy.astro.maneuver.transfer as transfer
import src.valladopy.constants as const

from ...conftest import custom_isclose


@pytest.fixture
def rinit_rfinal():
    rinit = const.RE + 191.3411
    rfinal = const.RE + 35781.34857
    return rinit, rfinal


@pytest.fixture
def nuinit_nufinal():
    return 0, np.pi


@pytest.mark.parametrize(
    "einit, efinal, deltava_expected, deltavb_expected",
    [
        # Vallado 2007, Example 6-1
        (0, 0, 2.457038635627854, 1.478187015430896),
        # Test non-zero eccentricities
        (0.00123, 0.024, 2.4522496326398677, 1.4410649593617622),
    ],
)
def test_hohmann(
    rinit_rfinal, nuinit_nufinal, einit, efinal, deltava_expected, deltavb_expected
):
    # Calculate Hohmann transfer
    deltava, deltavb, dtsec = transfer.hohmann(
        *rinit_rfinal, einit, efinal, *nuinit_nufinal
    )

    # Expected results
    assert custom_isclose(deltava, deltava_expected)
    assert custom_isclose(deltavb, deltavb_expected)
    assert custom_isclose(dtsec / const.MIN2SEC, 315.4027569935991)


@pytest.mark.parametrize(
    "einit, efinal, deltava_expected, deltavb_expected, deltavc_expected",
    [
        # Vallado 2007, Example 6-2
        (0, 0, 3.1562341430432026, 0.6773580299907941, 0.07046593708633986),
        # Test non-zero eccentricities
        (0.00123, 0.024, 3.151445140055216, 0.6773580299907941, 0.08278725521142083),
    ],
)
def test_bielliptic(
    rinit_rfinal,
    nuinit_nufinal,
    einit,
    efinal,
    deltava_expected,
    deltavb_expected,
    deltavc_expected,
):
    # Input values
    rinit, _ = rinit_rfinal
    rb = const.RE + 503873
    rfinal = const.RE + 376310

    # Calculate bielliptic transfer
    deltava, deltavb, deltavc, dtsec = transfer.bielliptic(
        rinit, rb, rfinal, einit, efinal, *nuinit_nufinal
    )

    # Expected results
    assert custom_isclose(deltava, deltava_expected)
    assert custom_isclose(deltavb, deltavb_expected)
    assert custom_isclose(deltavc, deltavc_expected)
    assert custom_isclose(dtsec / const.HR2SEC, 593.9199987732167)


@pytest.mark.parametrize(
    "efinal, deltavb_expected",
    [
        # Vallado 2007, Example 6-3
        (0, 2.1239372795577847),
        # Test non-zero eccentricity
        (0.024, 2.0814417965636647),
    ],
)
def test_one_tangent(rinit_rfinal, nuinit_nufinal, efinal, deltavb_expected):
    # Input values
    nuinit, _ = nuinit_nufinal
    nutran = np.radians(160)

    # Calculate one-tangent transfer
    deltava, deltavb, dtsec, etran, atran, vtrana, vtranb = transfer.one_tangent(
        *rinit_rfinal, efinal, nuinit, nutran
    )

    # Expected results
    assert custom_isclose(deltava, 2.575396079597792)
    assert custom_isclose(deltavb, deltavb_expected)
    assert custom_isclose(dtsec / const.MIN2SEC, 207.44540183723626)
    assert custom_isclose(etran, 0.7705727464030682)
    assert custom_isclose(atran, 28634.25027761329)
    assert custom_isclose(vtrana, 10.364786573669267)
    assert custom_isclose(vtranb, 2.2335537984044187)


def test_one_tangent_bad_transfer(rinit_rfinal, nuinit_nufinal):
    # Input values
    efinal = 0
    nuinit, _ = nuinit_nufinal
    nutran = np.radians(0)

    # Calculate one-tangent transfer
    with pytest.raises(ValueError):
        transfer.one_tangent(*rinit_rfinal, efinal, nuinit, nutran)


@pytest.mark.parametrize(
    "vinit, fpa, deltav_expected",
    [
        # Vallado 2007, Example 6-4
        (5.892311, 0, 1.5382018364126486),
        # Test non-zero flight path angle, Example 6-4
        (5.993824, np.radians(-6.79), 1.5537274747634255),
    ],
)
def test_incl_only(vinit, fpa, deltav_expected):
    # Input values
    deltai = np.radians(15)

    # Calculate inclination-only change
    deltav = transfer.incl_only(deltai, vinit, fpa)

    # Expected results
    assert custom_isclose(deltav, deltav_expected)


@pytest.mark.parametrize(
    "ecc, fpa, ifinal_exp, deltav_exp, arglat_init_exp, arglat_final_exp",
    [
        # Vallado 2007, Example 6-5
        (
            0,
            0,
            np.radians(55),
            3.694195175425934,
            np.radians(103.36472752868785),
            np.radians(76.63527247131213),
        ),
        # Test non-zero eccentricity flight path angle
        (
            0.0123,
            np.radians(11.23),
            np.radians(63.6589989314141),
            3.8891364302301996,
            np.radians(90),
            np.radians(77.4574164527168),
        ),
    ],
)
def test_node_only(ecc, fpa, ifinal_exp, deltav_exp, arglat_init_exp, arglat_final_exp):
    # Input values
    iinit = np.radians(55)
    deltaraan = np.radians(45)
    vinit = 5.892311
    incl = np.radians(55)

    # Calculate node change
    ifinal, deltav, arglat_init, arglat_final = transfer.node_only(
        iinit, ecc, deltaraan, vinit, fpa, incl
    )

    # Expected results
    assert custom_isclose(ifinal, float(ifinal_exp))
    assert custom_isclose(deltav, deltav_exp)
    assert custom_isclose(arglat_init, float(arglat_init_exp))
    assert custom_isclose(arglat_final, float(arglat_final_exp))


@pytest.mark.parametrize(
    "fpa, deltav_expected",
    [
        # Vallado 2007, Example 6-6
        (0, 3.615924548496319),
        # Test non-zero flight path angle
        (np.radians(11.23), 3.546691598327748),
    ],
)
def test_incl_and_node(fpa, deltav_expected):
    # Input values
    iinit = np.radians(55)
    ifinal = np.radians(40)
    deltaraan = np.radians(45)
    vinit = 5.892311

    # Calculate inclination and node change
    deltav, arglat_init, arglat_final = transfer.incl_and_node(
        iinit, ifinal, deltaraan, vinit, fpa
    )

    # Expected results
    assert custom_isclose(deltav, deltav_expected)
    assert custom_isclose(arglat_init, float(np.radians(128.9041397405515)))
    assert custom_isclose(arglat_final, float(np.radians(97.38034532660193)))


@pytest.mark.parametrize(
    "einit, efinal, deltai_init_exp, deltai_final_exp, deltava_exp, deltavb_exp, "
    "use_optimal",
    [
        # Vallado 2007, Table 6-3, optimal
        (
            0,
            0,
            np.radians(0.9173377342379527),
            np.radians(9.082662265762048),
            2.431758178760123,
            1.5091409647339542,
            True,
        ),
        # Vallado 2007, Table 6-3, non-optimal
        (
            0,
            0,
            np.radians(0.589638589451664),
            np.radians(9.410361410548337),
            2.4293288846983843,
            1.512146920521761,
            False,
        ),
        # Test non-zero eccentricities
        (
            0.00123,
            0.024,
            np.radians(0.9256322124753574),
            np.radians(9.074367787524643),
            2.4270919461253353,
            1.4724866802427696,
            True,
        ),
    ],
)
def test_min_combined(
    nuinit_nufinal,
    einit,
    efinal,
    deltai_init_exp,
    deltai_final_exp,
    deltava_exp,
    deltavb_exp,
    use_optimal,
):
    # Input values
    rinit = 6671.53  # km
    rfinal = 42163.95  # km
    iinit = np.radians(55)
    ifinal = np.radians(65)

    # Calculate minimum combined transfer
    deltai_init, deltai_final, deltava, deltavb, dtsec = transfer.min_combined(
        rinit,
        rfinal,
        einit,
        efinal,
        *nuinit_nufinal,
        iinit,
        ifinal,
        use_optimal=use_optimal,
    )

    # Expected results
    assert custom_isclose(deltai_init, float(deltai_init_exp))
    assert custom_isclose(deltai_final, float(deltai_final_exp))
    assert custom_isclose(deltava, deltava_exp)
    assert custom_isclose(deltavb, deltavb_exp)
    assert custom_isclose(dtsec / const.MIN2SEC, 316.4374908669237)


@pytest.mark.parametrize(
    "einit, nuinit, deltava_exp, gama_exp",
    [
        # Vallado 2007, Example 6-7
        (0, 0, 2.4696696299347445, np.radians(6.6314687256750915)),
        # Test non-zero eccentricity and true anomaly
        (0.00123, np.pi, 2.474429703935413, np.radians(6.618654628381083)),
    ],
)
def test_combined(einit, nuinit, deltava_exp, gama_exp):
    rinit = const.RE + 191
    rfinal = const.RE + 35780
    deltai = np.radians(-28.5)

    # Calculate combined transfer
    deltai1, deltai2, deltava, deltavb, dtsec, gama, gamb = transfer.combined(
        rinit, rfinal, einit, nuinit, deltai
    )

    # Expected results
    assert custom_isclose(deltai1, float(np.radians(-1.594961534092242)))
    assert custom_isclose(deltai2, float(np.radians(-26.90503846590776)))
    assert custom_isclose(deltava, deltava_exp)
    assert custom_isclose(deltavb, 1.8022162930889054)
    assert custom_isclose(dtsec / const.MIN2SEC, 315.3863523156078)
    assert custom_isclose(gama, float(gama_exp))
    assert custom_isclose(gamb, float(np.radians(50.53946267721802)))


@pytest.mark.parametrize(
    "rcstgt, einit, efinal, nufinal, kint, ktgt, phasef_exp, waittime_exp, deltav_exp",
    [
        # Vallado 2007, Example 6-8, Part 1
        (
            12756.274,
            0,
            0,
            0,
            1,
            1,
            np.radians(-20.0),
            225.6949830406384,
            -0.2192686615890107,
        ),
        # Vallado 2007, Example 6-8, Part 2
        (
            12756.274,
            0,
            0,
            0,
            2,
            2,
            np.radians(-20.0),
            464.6661415542555,
            -0.10648139790646738,
        ),
        # Vallado 2007, Example 6-9
        (
            42159.48,
            0,
            0,
            0,
            1,
            1,
            np.radians(85.391519024022),
            370.61383589095425,
            2.315710041662876,
        ),
        # Test non-zero eccentricities and final true anomaly
        (
            42159.48,
            0.00123,
            0.024,
            np.pi,
            1,
            1,
            np.radians(85.391519024022),
            370.61383589095425,
            2.2751512282382222,
        ),
    ],
)
def test_rendezvous_coplanar(
    rcstgt, einit, efinal, nufinal, kint, ktgt, phasef_exp, waittime_exp, deltav_exp
):
    # Input values
    rcsint = 12756.274
    phasei = np.radians(-20)
    nuinit = 0

    # Calculate rendezvous transfer
    phasef, waittime, deltav = transfer.rendezvous_coplanar(
        rcsint, rcstgt, phasei, einit, efinal, nuinit, nufinal, kint, ktgt
    )

    # Expected results
    assert custom_isclose(phasef, float(phasef_exp))
    assert custom_isclose(waittime / const.MIN2SEC, waittime_exp)
    assert custom_isclose(deltav, deltav_exp)


def test_rendezvous_coplanar_intersects():
    # Input values
    rcsint = rcstgt = 12756.274
    phasei = np.radians(-20)
    einit = efinal = nuinit = nufinal = 0
    kint, ktgt = 10, 1

    # Calculate rendezvous transfer
    with pytest.raises(ValueError):
        transfer.rendezvous_coplanar(
            rcsint, rcstgt, phasei, einit, efinal, nuinit, nufinal, kint, ktgt
        )


def test_rendezvous_noncoplanar():
    # Vallado, Example 6-10
    phasei = np.pi - np.radians(15)
    aint = 7143.51
    atgt = 42159.4855
    kint, ktgt = 1, 0
    nodeint = np.radians(45)
    truelon = np.radians(200)
    deltai = np.radians(-28.5)

    # Calculate rendezvous transfer
    out = transfer.rendezvous_noncoplanar(
        phasei, aint, atgt, kint, ktgt, nodeint, truelon, deltai
    )
    ttrans, tphase, dvphase, dvtrans1, dvtrans2, aphase = out

    # Expected results
    assert custom_isclose(ttrans / const.MIN2SEC, 320.9923637303245)
    assert custom_isclose(tphase / const.MIN2SEC, 450.7328392802076)
    assert custom_isclose(dvphase, 2.0762730792979385)
    assert custom_isclose(dvtrans1, 0.2226079943473298)
    assert custom_isclose(dvtrans2, 1.8024492966461878)
    assert custom_isclose(aphase, 19473.30414898291)
