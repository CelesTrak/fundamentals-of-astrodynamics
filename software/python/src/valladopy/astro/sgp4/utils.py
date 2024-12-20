from dataclasses import dataclass
from enum import Enum

import numpy as np


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
        xke = 60.0 / np.sqrt(radiusearthkm**3 / mu)
        j2 = 0.001082616
        j3 = -0.00000253881
        j4 = -0.00000165597
    elif wgs_model == WGSModel.WGS_84:
        mu = 398600.5
        radiusearthkm = 6378.137
        xke = 60.0 / np.sqrt(radiusearthkm**3 / mu)
        j2 = 0.00108262998905
        j3 = -0.00000253215306
        j4 = -0.00000161098761
    else:
        raise ValueError(f"Unknown option: {wgs_model}")

    return GravitationalConstants(1 / xke, mu, radiusearthkm, xke, j2, j3, j4, j3 / j2)
