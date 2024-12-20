import pytest


@pytest.fixture
def epoch():
    return 20630.3321544402


@pytest.fixture
def oe_params():
    ecc = 0.6877146
    incl = 1.11977881347003
    node = 4.87072001413786
    argp = 4.62102273937204
    n = 0.00874854701963024
    m = 0.353005058520617
    return ecc, incl, node, argp, n, m
