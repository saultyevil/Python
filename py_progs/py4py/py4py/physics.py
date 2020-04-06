"""
Physics functions

Things like calculating Keplerian velocity, doppler shifts and the like.
"""
import numpy as np
from astropy import constants as apc


def keplerian_velocity(mass: float, radius: float) -> float:
    """
    Calculates Keplerian velocity at given radius

    Args:
        mass (float): Object mass in kg
        radius (float): Orbital radius in m

    Returns:
        float: Orbital velocity in m/s
    """
    return np.sqrt(apc.G.value * mass / radius)


def doppler_shift_wave(line: float, vel: float) -> float:
    """
    Converts passed line and velocity into red/blue-shifted wavelength

    Args:
        line (float): Line wavelength (any length unit)
        vel (float): Doppler shift velocity (m/s)

    Returns:
        float: Doppler shifted line wavelength (as above)
    """
    return line * apc.c.value / (apc.c.value - vel)


def doppler_shift_vel(line: float, wave: float) -> float:
    """
    Converts passed red/blue-shifted wave into velocity

    Args:
        line (float):   Base line wavelength (any length unit)
        wave (float):   Doppler shifted line wavelength (as above)

    Returns:
        float:          Speed of Doppler shift
    """
    return apc.c.value * (wave - line) / line
