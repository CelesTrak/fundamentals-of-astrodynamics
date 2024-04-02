# -----------------------------------------------------------------------------
# Name: constants.py
# Author: David Vallado
# Date: 2 Apr 2007
#
# Copyright (c) 2024
# For license information, see LICENSE file
# -----------------------------------------------------------------------------

import numpy as np

###############################################################################
# Mathematical Operations
###############################################################################

SMALL = 1e-10

# Conversions
FT2M = 0.3048
MILE2M = 1609.344
NM2M = 1852
MILE2FT = 5280
MILEPH2KMPH = 0.44704
NMPH2KMPH = 0.5144444

# Time
DAY2SEC = 86400
MIN2SEC = 60


###############################################################################
# Astrodynamic Operations
###############################################################################

# Physical constants

# EGM-08 (Earth) constants used here
RE = 6378.1363                      # km
FLAT = 1.0 / 298.257223563
EARTHROT = 7.292115e-5              # rad/s
MU = 398600.4415                    # km^3/s^2
MUM = 3.986004415e14                # m^3/s^2

# Derived constants from the base values

# Earth eccentricity
ECCEARTH = np.sqrt(2.0 * FLAT - FLAT**2)
ECCEARTHSQRD = ECCEARTH**2

# Earth radius
RENM = RE / NM2M
REFT = RE * 1e3 / FT2M

# Orbital period
TUSEC = np.sqrt(RE**3 / MU)
TUMIN = TUSEC / MIN2SEC
TUDAY = TUSEC / DAY2SEC
TUDAYSID = TUSEC / 86164.090524

# Earth rotational angular velocity
OMEGAARTHPTU = EARTHROT * TUSEC
OMEGAARTHPMIN = EARTHROT * MIN2SEC

# Orbital velocity
VELKPS = np.sqrt(MU / RE)
VELFPS = VELKPS * 1e3 / FT2M
VELPDMIN = VELKPS * MIN2SEC / RE
DEGSEC = (180.0 / np.pi) / TUSEC
RADPDAY = 2.0 * np.pi * 1.002737909350795

# Astronomical distances & measurements
SPEEDOFLIGHT = 299792.458           # km/s
AU = 149597870.7                    # km
EARTH2MOON = 384400.0               # km
MOONRADIUS = 1738.0                 # km
SUNRADIUS = 696000.0                # km

# Masses in kg
MASSSUN = 1.9891e30
MASSEARTH = 5.9742e24
MASSMOON = 7.3483e22

# Standard gravitational parameters in km^3/s^2
MUSUN = 1.32712428e11
MUMOON = 4902.799
