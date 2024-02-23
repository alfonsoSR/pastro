from .base import Force
from ..types import CartesianState, CartesianStateDerivative


class ThirdBody(Force):

    def __init__(self, s_body: CartesianState, mu_body: float) -> None:

        self.body = s_body
        self.mu = mu_body

        return None

    def __call__(self, t: float, s: CartesianState,
                 fr: float = 0.) -> CartesianStateDerivative:

        # Position vector of body with respect to satellite
        s_sb = self.body - s

        # Parameters
        r_sb_mag = s_sb.r_mag
        r_sb_mag3 = r_sb_mag * r_sb_mag * r_sb_mag

        r_eb_mag = self.body.r_mag
        r_eb_mag3 = r_eb_mag * r_eb_mag * r_eb_mag

        p1 = self.mu / r_sb_mag3
        p2 = - self.mu / r_eb_mag3

        # Third body acceleration
        return CartesianStateDerivative(s.dx, s.dy, s.dz,
                                        p1 * s_sb.x + p2 * self.body.x,
                                        p1 * s_sb.y + p2 * self.body.y,
                                        p1 * s_sb.z + p2 * self.body.z)
