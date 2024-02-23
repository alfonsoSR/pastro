from .base import Force
from ..types import CartesianState, CartesianStateDerivative
import numpy as np


class SRP(Force):
    """Acceleration due to solar radiation pressure

    :param s_sun: Cartesian state vector of the Sun w.r.t Earth in equatorial
    inertial reference frame [m].
    :param Re: Earth radius [m]
    :param Rs: Sun radius [m] (Defaults to 6.96e8)
    :param Br: Radiation ballistic coefficient [m^2/kg] (Defaults to 0.02)
    :param Wsun: Solar constant [W/m^2] (Defaults to 1361)
    """

    def __init__(self, s_sun: CartesianState, Re: float, Rs: float = 6.96e8,
                 Br: float = 0.02, Wsun: float = 1361) -> None:

        self.sun = s_sun
        self.Re = Re
        self.Rs = Rs
        self.Br = Br
        self.c = 2.99792458e8               # Speed of light [m/s]
        self.P = Wsun / self.c

        return None

    def __call__(self, t: float, s: CartesianState,
                 fr: float = 0.) -> CartesianStateDerivative:

        # Position vector of satellite w.r.t Earth
        r_es_vec = s.r

        # Position vector of Sun w.r.t. Earth
        r_eo_vec = self.sun.r

        # Position vector of satellite w.r.t Sun
        r_os_vec = r_es_vec - r_eo_vec
        r_os_mag = np.linalg.norm(r_os_vec, axis=0)
        r_os_uvec = r_os_vec / r_os_mag

        # Acceleration not taking shadows into account
        common = self.P * self.Br * r_os_uvec

        if r_os_mag < self.sun.r_mag:
            return CartesianStateDerivative(s.dx, s.dy, s.dz,
                                            common[0], common[1], common[2])

        # Position vector of satellite w.r.t point P
        r_ps_vec = np.sum(r_os_uvec * r_es_vec, axis=0) * r_os_uvec
        r_ps_mag = np.linalg.norm(r_ps_vec, axis=0)

        # Position vector of point P w.r.t Earth
        r_ep_vec = r_es_vec - r_ps_vec
        r_ep_mag = np.linalg.norm(r_ep_vec, axis=0)

        # Height of point P above Earth's surface
        hg = r_ep_mag - self.Re

        # Apparent radius of the Sun at point P
        Rp = r_ps_mag * self.Rs / r_os_mag

        # Shadow function
        eta = hg / Rp

        if eta > 1.:
            return CartesianStateDerivative(s.dx, s.dy, s.dz,
                                            common[0], common[1], common[2])
        elif eta < -1.:
            return CartesianStateDerivative(s.dx, s.dy, s.dz, 0., 0., 0.)
        else:
            fs = (1. - (np.arccos(eta) / np.pi) +
                  (eta * np.sqrt(1. - eta * eta) / np.pi))
            return CartesianStateDerivative(s.dx, s.dy, s.dz, fs * common[0],
                                            fs * common[1], fs * common[2])
