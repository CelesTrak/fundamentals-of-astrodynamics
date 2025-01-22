import numpy as np

from ... import constants as const


# TODO: Move to a more general location


def time_of_flight(a: float) -> float:
    """Calculate the time of flight for a given semi-major axis.

    Args:
        a (float): Semi-major axis of the orbit in km

    Returns:
        float: Time of flight in seconds
    """
    return np.pi * np.sqrt(a**3 / const.MU)


def specific_mech_energy(a: float) -> float:
    """Computes specific mechanical energy at a given semi-major axis.

    Args:
        a (float): Semi-major axis of the orbit in km

    Returns:
        float: Specific mechanical energy
    """
    return -const.MU / (2 * a)


def velocity_mag(r: float, a: float) -> float:
    """Computes velocity magnitude at a given radius and specific mechanical energy.

    Args:
        r (float): Radius of the orbit in km
        a (float): Semi-major axis of the orbit in km

    Returns:
        float: Velocity magnitude in km/s
    """
    return np.sqrt(2 * ((const.MU / r) + specific_mech_energy(a)))


def angular_velocity(a: float) -> float:
    """Computes the angular velocity at a given semi-major axis.

    Args:
        a (float): Semi-major axis of the orbit in km

    Returns:
        float: Angular velocity in rad/s
    """
    return np.sqrt(const.MU / a**3)


def deltav(v1: float, v2: float, theta: float) -> float:
    """Computes the delta-v for a given two velocity magnitudes and an angle.

    Args:
        v1 (float): Initial velocity in km/s
        v2 (float): Final velocity in km/s
        theta (float): Angle in radians

    Returns:
        float: Delta-v in km/s
    """
    return np.sqrt(v1**2 + v2**2 - 2 * v1 * v2 * np.cos(theta))


def semimajor_axis(r: float, e: float, nu: float) -> float:
    """Computes the semi-major axis for a given radius, eccentricity, and true anomaly.

    Args:
        r (float): Radius of the orbit in km
        e (float): Eccentricity of the orbit
        nu (float): True anomaly in radians

    Returns:
        float: Semi-major axis in km
    """
    return (r * (1 + e * np.cos(nu))) / (1 - e**2)


def period(a: float) -> float:
    """Computes the period of an orbit given the semi-major axis.

    Args:
        a (float): Semi-major axis of the orbit in km

    Returns:
        float: Period of the orbit in seconds
    """
    return const.TWOPI * np.sqrt(a**3 / const.MU)


def lowuz(z: float) -> float:
    """Computes the control parameter u using Chebyshev polynomial approximation.

    References:
        Alfano: "Optimal Many-Revolution Orbit Transfer", Journal of Guidance, Vol 8,
                 No. 1, 1985, pp. 155-157.

    Args:
        z (float): Input parameter

    Returns:
        u (float): Computed control parameter (bounded between 0.000001 and 0.999999)
    """
    # Assign coefficients and starting values
    # fmt: off
    alpha = np.array([
        2.467410607, -1.907470562, 35.892442177, -214.672979624, 947.773273608,
        -2114.861134906, 2271.240058672, -1127.457440108, 192.953875268, 8.577733773
    ])

    beta = np.array([
        0.4609698838, 13.7756315324, -69.1245316678, 279.0671832500, -397.6628952136,
        -70.0139935047, 528.0334266841, -324.9303836520, 20.5838245170, 18.8165370778
    ])
    # fmt: on

    # Compute u using Chebyshev polynomial approximation
    alphasum, betasum, zterm = 0, 1, 1
    for i in range(10):
        zterm *= z
        alphasum += zterm * alpha[i]
        betasum += zterm * beta[i]

    u = abs(alphasum / betasum)

    # Clamp the value of u to be within the range [0.000001, 0.999999]
    u = max(0.000001, min(u, 0.999999))

    return u
