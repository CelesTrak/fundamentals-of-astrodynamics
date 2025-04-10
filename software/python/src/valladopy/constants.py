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

# Distances
KM2M = 1e3
FT2M = 0.3048
MILE2M = 1609.344
NM2M = 1852
MILE2FT = 5280
MILEPH2KMPH = 0.44704
NMPH2KMPH = 0.5144444

# Time
DAY2SEC = 86400
DAY2MIN = 1440
DAY2HR = 24
HR2SEC = 3600
MIN2SEC = 60
YR2DAY = 365.25
CENT2YR = 100
CENT2DAY = CENT2YR * YR2DAY

# Angles
HALFPI = np.pi / 2
TWOPI = 2 * np.pi
DEG2MIN = 60
DEG2ARCSEC = DEG2MIN * MIN2SEC
ARCSEC2RAD = np.radians(1 / DEG2ARCSEC)
DEG2SEC = np.degrees(TWOPI) / DAY2SEC
DEG2HR = np.degrees(TWOPI) / DAY2HR
HR2RAD = DEG2HR * np.radians(1)

###############################################################################
# Astrodynamic Operations
###############################################################################

# Time
J2000 = 2451545  # Julian date of the epoch J2000.0 (noon)
J2000_UTC = 2451544.5  # Julian date of the epoch J2000.0 in UTC (midnight)
JD_TO_MJD_OFFSET = 2400000.5  # offset between Julian and Modified Julian dates

# EGM-08 (Earth) constants used here
# fmt: off
RE = 6378.1363                      # km
FLAT = 1 / 298.257223563
EARTHROT = 7.292115e-5              # rad/s
MU = 398600.4415                    # km^3/s^2
MUM = 3.986004415e14                # m^3/s^2
J2 = 0.001082626174
J4 = -1.6198976e-06
# fmt: on

# Derived constants from the base values

# Sidereal day in seconds
SIDERALDAY_SEC = 86164.090524  # seconds

# Approximate Earth rotation
EARTHROT_APPROX = TWOPI / DAY2SEC  # rad/s

# Earth eccentricity
ECCEARTH = np.sqrt(2 * FLAT - FLAT**2)
ECCEARTHSQRD = ECCEARTH**2

# Earth radius
RENM = RE / NM2M
REFT = RE * 1e3 / FT2M

# Orbital period
TUSEC = np.sqrt(RE**3 / MU)
TUMIN = TUSEC / MIN2SEC
TUDAY = TUSEC / DAY2SEC
TUDAYSID = TUSEC / SIDERALDAY_SEC

# Earth rotation & rotational angular velocity
OMEGAARTHPTU = EARTHROT * TUSEC
OMEGAARTHPMIN = EARTHROT * MIN2SEC

# Orbital velocity
VELKPS = np.sqrt(MU / RE)
VELFPS = VELKPS * 1e3 / FT2M
VELPDMIN = VELKPS * MIN2SEC / RE
DEGSEC = (180 / np.pi) / TUSEC
RADPDAY = TWOPI * 1.002737909350795

# Astronomical distances & measurements
# fmt: off
SPEEDOFLIGHT = 299792.458           # km/s
AU2KM = 149597870.7                 # km
EARTH2MOON = 384400                 # km
MOONRADIUS = 1738                   # km
SUNRADIUS = 696000                  # km
# fmt: on

# Masses in kg
MASSSUN = 1.9891e30
MASSEARTH = 5.9742e24
MASSMOON = 7.3483e22

# Standard gravitational parameters in km^3/s^2
MUSUN = 1.32712428e11
MUMOON = 4902.799

# Obliquities
OBLIQUITYEARTH = np.radians(23.439291)
