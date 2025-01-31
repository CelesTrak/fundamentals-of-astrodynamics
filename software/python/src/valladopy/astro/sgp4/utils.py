# --------------------------------------------------------------------------------------
# Authors: David Vallado
# Date: 21 July 2006
#
# Copyright (c) 2024
# For license information, see LICENSE file
# --------------------------------------------------------------------------------------

from dataclasses import dataclass
from enum import Enum

import numpy as np
from pydantic import BaseModel

from ... import constants as const


class Classification(Enum):
    # fmt: off
    """Classification of the satellite."""
    Unclassified = "U"
    Classified = "C"


class SatRec(BaseModel):
    a: float = 0.0
    alta: float = 0.0
    altp: float = 0.0
    mo: float = 0.0
    mdot: float = 0.0
    argpo: float = 0.0
    argpdot: float = 0.0
    nodeo: float = 0.0
    nodedot: float = 0.0
    nodecf: float = 0.0
    cc1: float = 0.0
    cc4: float = 0.0
    cc5: float = 0.0
    omgcof: float = 0.0
    xmcof: float = 0.0
    eta: float = 0.0
    sinmao: float = 0.0
    delmo: float = 0.0
    d2: float = 0.0
    d3: float = 0.0
    d4: float = 0.0
    t2cof: float = 0.0
    t3cof: float = 0.0
    t4cof: float = 0.0
    t5cof: float = 0.0
    no: float = 0.0
    ecco: float = 0.0
    inclo: float = 0.0
    isimp: bool = False
    bstar: float = 0.0
    xfact: float = 0.0
    xlamo: float = 0.0
    atime: float = 0.0
    error: Enum = None
    t: float = 0.0
    aycof: float = 0.0
    xlcof: float = 0.0
    x1mth2: float = 0.0
    x7thm1: float = 0.0
    satnum: int = 0
    intldesg: str = ""
    epochyr: int = 0
    epochdays: float = 0.0
    ndot: float = 0.0
    nddot: float = 0.0
    elnum: int = 0
    revnum: int = 0
    no_kozai: float = 0.0
    jdsatepoch: float = 0.0
    jdsatepochf: float = 0.0
    init: bool = True
    classification: Classification = None


class WGSModel(Enum):
    WGS_72_LOW_PRECISION = 721
    WGS_72 = 72
    WGS_84 = 84


@dataclass
class GravitationalConstants:
    tumin: float
    mu: float
    radiusearthkm: float
    xke: float
    j2: float
    j3: float
    j4: float
    j3oj2: float


def getgravc(wgs_model: WGSModel) -> GravitationalConstants:
    """Returns the gravitational constants based on the specified WGS model.

    References:
        - NORAD SpaceTrack Report #3
        - Vallado, Crawford, Hujsak, Kelso, 2006

    Args:
        wgs_model (WGSModel): The WGS model to use

    Returns:
        GravitationalConstants: A data structure containing the gravitational constants
    """
    if wgs_model == WGSModel.WGS_72_LOW_PRECISION:
        mu = 398600.79964
        radiusearthkm = 6378.135
        xke = 0.0743669161
        j2 = 0.001082616
        j3 = -0.00000253881
        j4 = -0.00000165597
    elif wgs_model == WGSModel.WGS_72:
        mu = 398600.8
        radiusearthkm = 6378.135
        xke = const.MIN2SEC / np.sqrt(radiusearthkm**3 / mu)
        j2 = 0.001082616
        j3 = -0.00000253881
        j4 = -0.00000165597
    elif wgs_model == WGSModel.WGS_84:
        mu = 398600.5
        radiusearthkm = 6378.137
        xke = const.MIN2SEC / np.sqrt(radiusearthkm**3 / mu)
        j2 = 0.00108262998905
        j3 = -0.00000253215306
        j4 = -0.00000161098761
    else:
        raise ValueError(f"Unknown option: {wgs_model}")

    return GravitationalConstants(1 / xke, mu, radiusearthkm, xke, j2, j3, j4, j3 / j2)
