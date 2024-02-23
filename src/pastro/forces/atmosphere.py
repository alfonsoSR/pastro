import numpy as np


class DensityModel:
    """Base class for atmospheric density models"""

    def __init__(self, *args) -> None:
        raise NotImplementedError

    def __call__(self, h: float) -> float:
        raise NotImplementedError


class ExponentialAtmosphere(DensityModel):

    def __init__(self, rho0: float, Hs: float) -> None:

        self.rho0 = rho0
        self.beta = 1. / Hs

        return None

    def __call__(self, h: float) -> float:

        return self.rho0 * np.exp(- h * self.beta)
