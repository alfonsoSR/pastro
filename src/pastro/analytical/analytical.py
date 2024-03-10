from ..types import (
    KeplerianState,
    Vector,
    time2mean_anomaly,
    true2eccentric_anomaly,
    eccentric2mean_anomaly,
    eccentric2true_anomaly,
)
import numpy as np
from ..root_finding import enrke


def keplerian_orbit(s0: KeplerianState, epochs: Vector, mu: float) -> KeplerianState:

    if not s0.is_scalar:
        raise ValueError("Initial state must be scalar")

    # Initial value of mean anomaly
    E0 = true2eccentric_anomaly(s0.nu, s0.e)
    M0 = eccentric2mean_anomaly(E0, s0.e)

    # True anomaly as function of time
    assert isinstance(M0, float)
    M = time2mean_anomaly(epochs, s0.a, mu, M0)
    E = enrke(M, s0.e)
    nu_wrapped = eccentric2true_anomaly(E, s0.e)
    assert isinstance(nu_wrapped, np.ndarray)
    nu = np.unwrap(nu_wrapped, period=360.0)
    assert np.allclose(nu[0], s0.nu)

    # Generate keplerian state
    base = np.ones_like(epochs, dtype=np.float64)
    return KeplerianState(
        s0.a * base,
        s0.e * base,
        s0.i * base,
        s0.Omega * base,
        s0.omega * base,
        nu,
    )


# def keplerian_orbit(s0: KeplerianState, epochs: Vector, mu: float) -> KeplerianState:
#     """Generate an ideal keplerian trajectory

#     :param s0: Keplerian elements at initial epoch
#     :param epochs: Epoch at which to calculate the state [s]
#     :param mu: Standard gravitational parameter of central body [m^3/s^2]
#     """

#     a = np.ones_like(epochs) * s0.a
#     e = np.ones_like(epochs) * s0.e
#     i = np.ones_like(epochs) * s0.i
#     Omega = np.ones_like(epochs) * s0.Omega
#     omega = np.ones_like(epochs) * s0.omega

#     # Calculate mean anomaly as function of time
#     # M = (t - t0) * sqrt(mu / a^3)
#     M = (epochs - epochs[0]) * np.sqrt(mu / (s0.a * s0.a * s0.a))

#     # Calculate eccentric anomaly from mean anomaly
#     # M = E - e * sin(E)
#     # Newton-Raphson using M as initial guess
#     def f(E, e, M):
#         return E - e * np.sin(E) - M

#     def fprime(E, e, M):
#         return 1 - e * np.cos(E)

#     def fprime2(E, e, M):
#         return e * np.sin(E)

#     guess_E = M

#     E = newton(f, guess_E, fprime, fprime2=fprime2, args=(e, M))

#     # Calculate true anomaly from eccentric anomaly
#     nu = np.rad2deg(2.0 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2)))

#     return KeplerianState(a, e, i, Omega, omega, nu)
