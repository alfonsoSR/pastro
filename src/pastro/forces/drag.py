from .base import Force
from ..types import CartesianState, CartesianStateDerivative
from .atmosphere import DensityModel
import numpy as np


class Drag(Force):
    """Acceleration due to atmospheric drag

    :param Bd: Ballistic coefficient [m^2/kg] (Defaults to 0.01)
    :param atmos: Atmospheric density model
    :param omega: Angular velocity of Earth [rad/s]
    """

    def __init__(self, atmos: DensityModel, omega: float,
                 Re: float, Bd: float = 0.01) -> None:

        self.Bd = Bd
        self.rho = atmos
        self.omega_vec = np.array([0., 0., omega])
        self.Re = Re

        return None

    def __call__(self, t: float, s: CartesianState,
                 fr: float = 0.) -> CartesianStateDerivative:

        # Height over surface of Earth
        h = s.r_mag - self.Re
        assert (isinstance(h, float))

        # Velocity with respect to rotating atmosphere
        u_vec = s.u - np.cross(self.omega_vec, s.r, axis=0)
        u_mag = np.linalg.norm(u_vec, axis=0)

        # Acceleration due to atmospheric drag
        a_drag = -0.5 * self.rho(h) * self.Bd * u_mag * u_vec

        return CartesianStateDerivative(s.dx, s.dy, s.dz,
                                        a_drag[0], a_drag[1], a_drag[2])
