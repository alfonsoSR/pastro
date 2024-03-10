from .core import PlotSetup, GenericPlot, MultiPlot
from ..types import Vector, CartesianState, KeplerianState
from .core import PlotSetup, SinglePlot, MultiPlot

# TODO: FIX OVERRIDE ISSUE IN PLOT METHODS


class PlotKeplerianState(MultiPlot):
    """Plot evolution of classical keplerian elements"""

    def __init__(self, setup: PlotSetup) -> None:

        if setup.xlabel is None:
            _xlabel = "Days past initial epoch"
        else:
            _xlabel = setup.xlabel

        setup.subplots = (3, 2)

        self.a_setup = PlotSetup(ylabel="$a$ [km]", xlabel=_xlabel)
        self.e_setup = PlotSetup(ylabel="$e$", xlabel=_xlabel)
        self.i_setup = PlotSetup(ylabel="$i$ [deg]", xlabel=_xlabel)
        self.omega_setup = PlotSetup(ylabel="AoP [deg]", xlabel=_xlabel)
        self.Omega_setup = PlotSetup(ylabel="RAAN [deg]", xlabel=_xlabel)
        self.nu_setup = PlotSetup(ylabel="TA [deg]", xlabel=_xlabel)

        super().__init__(setup)

        return None

    def plot(self, time: Vector, state: KeplerianState) -> str:
        """Plot evolution of classical keplerian elements

        :param time: Time vector [s]
        :param state: Keplerian state
        """

        # Turn time vector into days past initial epoch
        dt = (time - time[0]) / (24.0 * 3600.0)

        self.add_plot(self.a_setup).plot(dt, state.a * 1e-3)
        self.add_plot(self.e_setup).plot(dt, state.e)
        self.add_plot(self.i_setup).plot(dt, state.i)
        self.add_plot(self.omega_setup).plot(dt, state.omega)
        self.add_plot(self.Omega_setup).plot(dt, state.Omega)
        self.add_plot(self.nu_setup).plot(dt, state.nu)

        return self.__call__()


class PlotCartesianState(MultiPlot):

    def __init__(self, setup: PlotSetup) -> None:

        if setup.xlabel is None:
            _xlabel = "Days past initial epoch"
        else:
            _xlabel = setup.xlabel

        setup.subplots = (3, 2)

        self.x_setup = PlotSetup(ylabel="$x$ [km]", xlabel=_xlabel)
        self.y_setup = PlotSetup(ylabel="$y$ [km]", xlabel=_xlabel)
        self.z_setup = PlotSetup(ylabel="$z$ [km]", xlabel=_xlabel)
        self.vx_setup = PlotSetup(ylabel="$v_x$ [km/s]", xlabel=_xlabel)
        self.vy_setup = PlotSetup(ylabel="$v_y$ [km/s]", xlabel=_xlabel)
        self.vz_setup = PlotSetup(ylabel="$v_z$ [km/s]", xlabel=_xlabel)

        super().__init__(setup)

        return None

    def plot(self, time: Vector, state: CartesianState) -> str:
        """Plot evolution of classical keplerian elements

        :param time: Time vector [s]
        :param state: Keplerian state
        """

        # Turn time vector into days past initial epoch
        dt = (time - time[0]) / (24.0 * 3600.0)

        self.add_plot(self.x_setup).plot(dt, state.x * 1e-3)
        self.add_plot(self.y_setup).plot(dt, state.y * 1e-3)
        self.add_plot(self.z_setup).plot(dt, state.z * 1e-3)
        self.add_plot(self.vx_setup).plot(dt, state.dx * 1e-3)
        self.add_plot(self.vy_setup).plot(dt, state.dy * 1e-3)
        self.add_plot(self.vz_setup).plot(dt, state.dz * 1e-3)

        return self.__call__()


class CompareCartesianOrbits(MultiPlot):

    def __init__(self, setup: PlotSetup) -> None:

        if setup.xlabel is not None:
            _xlabel = setup.xlabel
        else:
            _xlabel = "Days past initial epoch"

        setup.subplots = (3, 2)

        self.x_setup = PlotSetup(ylabel=r"$\Delta x$ [km]", xlabel=_xlabel)
        self.y_setup = PlotSetup(ylabel=r"$\Delta y$ [km]", xlabel=_xlabel)
        self.z_setup = PlotSetup(ylabel=r"$\Delta z$ [km]", xlabel=_xlabel)
        self.vx_setup = PlotSetup(ylabel=r"$\Delta \dot x$ [km/s]", xlabel=_xlabel)
        self.vy_setup = PlotSetup(ylabel=r"$\Delta \dot y$ [km/s]", xlabel=_xlabel)
        self.vz_setup = PlotSetup(ylabel=r"$\Delta \dot z$ [km/s]", xlabel=_xlabel)

        super().__init__(setup)

        return None

    def plot(
        self, time: Vector, orbit: CartesianState, reference: CartesianState
    ) -> str:
        """Plot difference between two sets of cartesian state vectors

        :param time: Time vector [s]
        :param orbit: Cartesian states of orbit.
        :param reference: Cartesian states of reference orbit.
        """

        # Turn time vector into days past initial epoch
        dt = (time - time[0]) / (24.0 * 3600.0)
        ds = orbit - reference

        self.add_plot(self.x_setup).plot(dt, ds.x * 1e-3)
        self.add_plot(self.y_setup).plot(dt, ds.y * 1e-3)
        self.add_plot(self.z_setup).plot(dt, ds.z * 1e-3)
        self.add_plot(self.vx_setup).plot(dt, ds.dx * 1e-3)
        self.add_plot(self.vy_setup).plot(dt, ds.dy * 1e-3)
        self.add_plot(self.vz_setup).plot(dt, ds.dz * 1e-3)

        return self.__call__()


class CompareKeplerianOrbits(MultiPlot):

    def __init__(self, setup: PlotSetup) -> None:

        if setup.xlabel is not None:
            _xlabel = setup.xlabel
        else:
            _xlabel = "Days past initial epoch"

        setup.subplots = (3, 2)

        self.a_setup = PlotSetup(ylabel=r"$\Delta a$ [km]", xlabel=_xlabel)
        self.e_setup = PlotSetup(ylabel=r"$\Delta e$", xlabel=_xlabel)
        self.i_setup = PlotSetup(ylabel=r"$\Delta i$ [deg]", xlabel=_xlabel)
        self.omega_setup = PlotSetup(ylabel=r"$\Delta \omega$ [deg]", xlabel=_xlabel)
        self.Omega_setup = PlotSetup(ylabel=r"$\Delta \Omega$ [deg]", xlabel=_xlabel)
        self.nu_setup = PlotSetup(ylabel=r"$\Delta \nu$ [deg]", xlabel=_xlabel)

        super().__init__(setup)

        return None

    def plot(
        self, time: Vector, orbit: KeplerianState, reference: KeplerianState
    ) -> str:
        """Plot difference between two sets of keplerian state vectors

        :param time: Time vector [s]
        :param orbit: Keplerian states of orbit.
        :param reference: Keplerian states of reference orbit.
        """

        # Turn time vector into days past initial epoch
        dt = (time - time[0]) / (24.0 * 3600.0)
        ds = orbit - reference

        self.add_plot(self.a_setup).plot(dt, ds.a * 1e-3)
        self.add_plot(self.e_setup).plot(dt, ds.e)
        self.add_plot(self.i_setup).plot(dt, ds.i)
        self.add_plot(self.omega_setup).plot(dt, ds.omega)
        self.add_plot(self.Omega_setup).plot(dt, ds.Omega)
        self.add_plot(self.nu_setup).plot(dt, ds.nu)

        return self.__call__()


# class PlotCartesianState(MultiPlot):
#     """Plot evolution of cartesian state components"""

#     def __init__(self, setup: PlotSetup) -> None:

#         if setup.xlabel is None:
#             _xlabel = "Days past initial epoch"
#         else:
#             _xlabel = setup.xlabel

#         setup.subplots = (3, 2)

#         x_setup = PlotSetup(ylabel="$x$ [km]", xlabel=_xlabel)
#         y_setup = PlotSetup(ylabel="$y$ [km]", xlabel=_xlabel)
#         z_setup = PlotSetup(ylabel="$z$ [km]", xlabel=_xlabel)
#         vx_setup = PlotSetup(ylabel="$v_x$ [km/s]", xlabel=_xlabel)
#         vy_setup = PlotSetup(ylabel="$v_y$ [km/s]", xlabel=_xlabel)
#         vz_setup = PlotSetup(ylabel="$v_z$ [km/s]", xlabel=_xlabel)

#         subplot_setups: list[PlotSetup] = [
#             x_setup,
#             y_setup,
#             z_setup,
#             vx_setup,
#             vy_setup,
#             vz_setup,
#         ]

#         super().__init__(setup, subplot_setups, GenericPlot)

#         return None

#     def plot(self, time: Vector, state: CartesianState) -> None:
#         """Plot evolution of classical keplerian elements

#         :param time: Time vector [s]
#         :param state: Keplerian state
#         """

#         # Turn time vector into days past initial epoch
#         dt = (time - time[0]) / (24.0 * 3600.0)

#         self.generators[0].plot(dt, state.x * 1e-3)  # x [m] -> [km]
#         self.generators[1].plot(dt, state.y * 1e-3)  # y [m] -> [km]
#         self.generators[2].plot(dt, state.z * 1e-3)  # z [m] -> [km]
#         self.generators[3].plot(dt, state.dx * 1e-3)  # dx [m/s] -> [km/s]
#         self.generators[4].plot(dt, state.dy * 1e-3)  # dy [m/s] -> [km/s]
#         self.generators[5].plot(dt, state.dz * 1e-3)  # dz [m/s] -> [km/s]

#         self._post()

#         return None


# class CompareCartesian(MultiPlot):
#     """Compare two sets of cartesian states"""

#     def __init__(self, setup: PlotSetup) -> None:

#         if setup.xlabel is None:
#             _xlabel = "Days past initial epoch"
#         else:
#             _xlabel = setup.xlabel

#         setup.subplots = (3, 2)

#         x_setup = PlotSetup(ylabel=r"$\Delta x$ [km]", xlabel=_xlabel)
#         y_setup = PlotSetup(ylabel=r"$\Delta y$ [km]", xlabel=_xlabel)
#         z_setup = PlotSetup(ylabel=r"$\Delta z$ [km]", xlabel=_xlabel)
#         dx_setup = PlotSetup(ylabel=r"$\Delta \dot x$ [km/s]", xlabel=_xlabel)
#         dy_setup = PlotSetup(ylabel=r"$\Delta \dot y$ [km/s]", xlabel=_xlabel)
#         dz_setup = PlotSetup(ylabel=r"$\Delta \dot z$ [km/s]", xlabel=_xlabel)

#         subplot_setups: list[PlotSetup] = [
#             x_setup,
#             y_setup,
#             z_setup,
#             dx_setup,
#             dy_setup,
#             dz_setup,
#         ]

#         super().__init__(setup, subplot_setups, GenericPlot)

#         return None

#     def plot(
#         self, time: Vector, target: CartesianState, reference: CartesianState
#     ) -> None:
#         """Plot difference between two sets of cartesian state vectors

#         :param time: Time vector [s]
#         :param target: Cartesian states of orbit.
#         :param reference: Cartesian states of reference orbit.
#         """

#         # Turn time vector into days past initial epoch
#         dt = (time - time[0]) / (24.0 * 3600.0)

#         self.generators[0].plot(dt, (target.x - reference.x) * 1e-3)
#         self.generators[1].plot(dt, (target.y - reference.y) * 1e-3)
#         self.generators[2].plot(dt, (target.z - reference.z) * 1e-3)
#         self.generators[3].plot(dt, (target.dx - reference.dx) * 1e-3)
#         self.generators[4].plot(dt, (target.dy - reference.dy) * 1e-3)
#         self.generators[5].plot(dt, (target.dz - reference.dz) * 1e-3)

#         self._post()

#         return None


class PlotRVMagnitudes(MultiPlot):
    """Plot the magnitude of the position and velocity vectors"""

    def __init__(self, setup: PlotSetup) -> None:

        if setup.xlabel is None:
            _xlabel = "Days past initial epoch"
        else:
            _xlabel = setup.xlabel

        setup.subplots = (2, 1)

        r_setup = PlotSetup(ylabel=r"$\| \vec{r} \|$ [km]", xlabel=_xlabel)
        v_setup = PlotSetup(ylabel=r"$\| \vec{v} \|$ [km/s]", xlabel=_xlabel)

        subplot_setups: list[PlotSetup] = [r_setup, v_setup]
        super().__init__(setup, subplot_setups, GenericPlot)

        return None

    def plot(self, time: Vector, state: CartesianState) -> None:
        """Plot magnitude of position and velocity vectors

        :param time: Time vector [s]
        :param state: Cartesian states of orbit.
        """

        # Turn time vector into days past initial epoch
        dt = (time - time[0]) / (24.0 * 3600.0)

        self.generators[0].plot(dt, state.r_mag * 1e-3)  # r [m] -> [km]
        self.generators[1].plot(dt, state.u_mag * 1e-3)  # v [m/s] -> [km/s]

        self._post()

        return None


# NOTE: NOT WORKING PROPERLY
# class PlotGroundTrack(MultiPlot):
#     """Plot ground track and orbital radius"""

#     def __init__(self, setup: PlotSetup) -> None:

#         if setup.xlabel is None:
#             _xlabel = "Days past initial epoch"
#         else:
#             _xlabel = setup.xlabel

#         setup.subplots = (2, 1)

#         r_setup = PlotSetup(ylabel="Orbital radius [km]", xlabel=_xlabel)
#         track_setup = PlotSetup(ylabel="Latitude [deg]",
#                                 xlabel="Longitude [deg]")

#         subplot_setups = [r_setup, track_setup]
#         super().__init__(setup, subplot_setups, SinglePlot)

#         return None

#     def plot(self, time: Vector, state: SphericalGeocentric):
#         """Plot ground track and orbital radius

#         :param time: Time vector [s]
#         :param state: Spherical geocentric states of orbit.
#         """

#         # Turn time vector into days past initial epoch
#         dt = (time - time[0]) / (24. * 3600.)

#         # Remove discontinuities in longitude
#         assert isinstance(state.lon, np.ndarray)
#         diff = (state.lon[:-1] - state.lon[1:]) > 190.
#         diff = np.concatenate((diff, [False]))
#         state.lon[diff] = np.nan

#         self.generators[0].plot(dt, state.r * 1e-3)    # r [m] -> [km]
#         self.generators[1].plot(state.lon, state.lat)
#         self.axes[1].set_xlim(-180, 181)

#         self._post()

#         return None
