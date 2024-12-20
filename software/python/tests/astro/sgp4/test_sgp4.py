import src.valladopy.astro.sgp4.sgp4 as sgp4

from ...conftest import custom_isclose


def test_sgp4(epoch, oe_params):
    # Inputs
    ecco, inclo, *_ = oe_params
    xke = 0.0743669161331734
    j2 = 0.001082616
    no_kozai = 0.00874808688806747

    # Expected outputs
    expected = {
        "ainv": 0.24008817584663986,
        "ao": 4.165136398215487,
        "con41": -0.43002188663163776,
        "con42": 0.05003647771939623,
        "cosio": 0.43588152571096744,
        "cosio2": 0.18999270445612076,
        "eccsq": 0.47295137105315993,
        "omeosq": 0.5270486289468401,
        "posq": 4.819032241803303,
        "rp": 1.3007112861712828,
        "rteosq": 0.7259811491676903,
        "sinio": 0.9000040530708066,
        "gsto": 0.574180126902192,
        "no_unkozai": 0.008748547019630244,
    }

    # Call method
    sgp4init_out = sgp4.initl(xke, j2, epoch, ecco, inclo, no_kozai)

    # Check results
    for key in expected:
        assert custom_isclose(getattr(sgp4init_out, key), expected[key])
