from .base import Force
from ..types import CartesianState, CartesianStateDerivative


class PointMass(Force):

    def __init__(self, mu: float) -> None:

        self.mu = mu

        return None

    def __call__(self, t: float, s: CartesianState,
                 fr: float = 0.) -> CartesianStateDerivative:

        f = - self.mu / (s.r_mag * s.r_mag * s.r_mag)

        return CartesianStateDerivative(s.dx, s.dy, s.dz,
                                        f * s.x, f * s.y, f * s.z)


class J2(Force):

    def __init__(self, mu: float, j2: float, R: float) -> None:

        self.mu = mu
        self.j2 = j2
        self.R = R

        return None

    def __call__(self, t: float, s: CartesianState,
                 fr: float = 0.) -> CartesianStateDerivative:

        r = s.r_mag
        r2 = r * r
        r5 = r2 * r2 * r

        p1 = -1.5 * self.mu * self.j2 * self.R * self.R / r5
        p2 = -5. * s.z * s.z / r2

        return CartesianStateDerivative(s.dx, s.dy, s.dz,
                                        p1 * s.x * (1. + p2),
                                        p1 * s.y * (1. + p2),
                                        p1 * s.z * (3. + p2))
