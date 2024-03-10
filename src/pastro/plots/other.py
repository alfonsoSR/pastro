from .core import SinglePlot, BasePlot, MultiPlot, GenericPlot, PlotSetup
from ..types import Vector, CartesianStateDerivative
import numpy as np


class PlotAcceleration(MultiPlot):

    def __init__(self, setup: PlotSetup = PlotSetup()) -> None:

        setup.subplots = (3, 1)

        super().__init__(setup)

        if setup.xlabel is None:
            self.xlabel = "Days past initial epoch"
        else:
            self.xlabel = setup.xlabel

        self.x_setup = PlotSetup(ylabel=r"$\ddot x$ [m/s^2]", xlabel=self.xlabel)
        self.y_setup = PlotSetup(ylabel=r"$\ddot y$ [m/s^2]", xlabel=self.xlabel)
        self.z_setup = PlotSetup(ylabel=r"$\ddot z$ [m/s^2]", xlabel=self.xlabel)

        return None

    def plot(
        self, epochs: Vector, acceleration: CartesianStateDerivative | Vector
    ) -> str:

        if isinstance(acceleration, np.ndarray):
            if acceleration.shape[0] != 3:
                raise ValueError("Acceleration must be a 3xN array")
            zero = np.zeros(acceleration.shape[1])
            ds = CartesianStateDerivative(
                zero, zero, zero, acceleration[0], acceleration[1], acceleration[2]
            )
        else:
            ds = acceleration

        dt = (epochs - epochs[0]) / (24.0 * 3600.0)

        self.add_plot(self.x_setup).plot(dt, ds.ddx)
        self.add_plot(self.y_setup).plot(dt, ds.ddy)
        self.add_plot(self.z_setup).plot(dt, ds.ddz)

        return self.__call__()
