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
