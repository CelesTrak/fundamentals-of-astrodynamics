import pytest

from src.valladopy.astro.sgp4.utils import SatRec

from ...conftest import custom_isclose


def test_dscom(ds, dscom_data):
    # Call method
    ds.dscom(tc=0)

    # Check results
    for key in dscom_data.__dataclass_fields__:
        assert custom_isclose(getattr(ds.dscom_out, key), getattr(dscom_data, key))


def test_dpper(ds, dscom_data):
    # Set dscom_out
    ds.dscom_out = dscom_data

    # Call method
    ds.dpper(t=0)

    # Check results
    assert custom_isclose(ds.ep, 0.6871280305587458)
    assert custom_isclose(ds.inclp, 1.1200831160846698)
    assert custom_isclose(ds.nodep, 4.869997828497486)
    assert custom_isclose(ds.argpp, 4.62188106182903)
    assert custom_isclose(ds.mp, 0.35183574006201346)


def test_dsinit(ds, dscom_data, dsinit_data):
    # Set dscom_out
    ds.dscom_out = dscom_data

    # Inputs
    satrec = SatRec(
        ecco=ds.ep,
        nodeo=ds.nodep,
        argpo=ds.argpp,
        no=ds.np_,
        mo=ds.mp,
        nodedot=-1.28456721580123e-06,
        mdot=0.00874808688663313,
    )
    tc = 0
    xke = 0.0743669161331734
    gsto = 0.574180126924752
    xpidot = -1.35893730845456e-06
    eccsq = dscom_data.emsq
    nodem = argpm = mm = 0
    inclm = ds.inclp

    # Call method
    ds.dsinit(satrec, xke, tc, gsto, xpidot, eccsq, inclm, nodem, argpm, mm)

    # Check results
    for key in dsinit_data.__dataclass_fields__:
        assert custom_isclose(getattr(ds.dsinit_out, key), getattr(dsinit_data, key))


@pytest.mark.parametrize(
    "irez, nm_expected, mm_expected, dndt_expected",
    [
        (0, 0.00874854701963024, 1.4027866098606812, 0),
        (1, 0.00874854701963024, -5.730064714256859, 0),
        (2, 0.008748547084955744, -4.880398693573535, 6.532550368698598e-11),
    ],
)
def test_dspace(dsinit_data, ds, irez, nm_expected, mm_expected, dndt_expected):
    # Inputs
    tc, gsto = 120, 0.574180126924752
    satrec = SatRec(t=tc, argpo=ds.argpp, argpdot=-7.43700926533354e-08, no=ds.np_)

    # Update DeepSpace object with more inputs
    ds.dsinit_out = dsinit_data
    ds.dsinit_out.irez = irez
    ds.dsinit_out.argpm = 4.62101381496092
    ds.dsinit_out.nodem = 4.87056586607196
    ds.dsinit_out.mm = 1.40277548491659

    # Call method
    ds.dspace(satrec, tc, gsto)

    # Check results
    assert custom_isclose(ds.dsinit_out.atime, 0)
    assert custom_isclose(ds.dsinit_out.em, 0.6877111337428518)
    assert custom_isclose(ds.dsinit_out.inclm, 1.1197776393770376)
    assert custom_isclose(ds.dsinit_out.argpm, 4.621012116335631)
    assert custom_isclose(ds.dsinit_out.nodem, 4.870558516883772)
    assert custom_isclose(ds.dsinit_out.nm, nm_expected)
    assert custom_isclose(ds.dsinit_out.mm, mm_expected)
    assert custom_isclose(ds.dsinit_out.dndt, dndt_expected)
    assert custom_isclose(ds.dsinit_out.xli, 2.66289952576725)
    assert custom_isclose(ds.dsinit_out.xni, 0.00874854701963024)
