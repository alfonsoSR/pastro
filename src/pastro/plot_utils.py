import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import axisartist
from mpl_toolkits.axes_grid1 import host_subplot
from pathlib import Path
from typing import Sequence, overload
from dataclasses import dataclass
from .utils import global_diff, cartesian2spherical, get_acceleration
from .types import (Vector, CartesianState,
                    KeplerianState, SphericalGeocentric)
from .forces import ForceModel


@dataclass
class PlotConfig:
    """Plot configuration.
    :param figsize: Size of the figure.
    :param title: Title of the figure.
    :param xlabel: Label for the x-axis.
    :param ylabel: Label for the y-axis.
    :param right_label: Label for the right y-axis.
    :param parasite_label: Label for the parasite y-axis.
    :param yscale: Scale for the y-axis.
    :param xscale: Scale for the x-axis.
    :param legend: If True, show legend.
    :param show: If True, show the plot.
    :param save: If True, save the plot.
    :param save_dir: Directory to save the plot.
    :param save_name: Name of the file to save.
    """

    figsize: tuple[int, int] = (10, 4)
    title: str | None = None
    xlabel: str | None = None
    ylabel: str | None = None
    right_label: str | None = None
    parasite_label: str | None = None
    yscale: str | None = None
    xscale: str | None = None
    legend: bool = False
    show: bool = True
    save: bool = False
    save_dir: str = "plots"
    save_name: str = ""


def __save_and_show(fig, config) -> str:

    path = ""
    if config.save:
        Path(config.save_dir).mkdir(parents=True, exist_ok=True)
        path = f"{config.save_dir}/{config.save_name}"
        fig.savefig(path)

    if config.show:
        plt.show()

    return path


class SimplePlot:
    """Base class for a figure with a single plot"""

    def __init__(self, config: PlotConfig) -> None:

        self.config = config

        self.fig = plt.figure(figsize=self.config.figsize)
        self.ax = host_subplot(
            111, axes_class=axisartist.Axes, figure=self.fig)    # type: ignore

        if self.config.title is not None:
            self.fig.suptitle(self.config.title)

        if self.config.xlabel is not None:
            self.ax.set_xlabel(self.config.xlabel)

        if self.config.ylabel is not None:
            self.ax.set_ylabel(self.config.ylabel)

        if self.config.yscale is not None:
            self.ax.set_yscale(self.config.yscale)

        if self.config.xscale is not None:
            self.ax.set_xscale(self.config.xscale)

        self.path = ""

        return None
        
    def post(self) -> None:
        
        if self.config.legend:
            self.ax.legend()
            
        self.fig.tight_layout()
        self.fig.subplots_adjust(left=0.1, right=0.9)

        if self.config.save:
            Path(self.config.save_dir).mkdir(parents=True, exist_ok=True)
            self.path = f"{self.config.save_dir}/{self.config.save_name}"
            self.fig.savefig(self.path)

        if self.config.show:
            plt.show(block=True)

        plt.close(self.fig)

        return None

    def default(self, x, y) -> str:

        self.ax.plot(x, y)
        self.post()

        return self.path


class DoubleAxis(SimplePlot):
    
    def __init__(self, config: PlotConfig) -> None:
        
        super().__init__(config)
        
        # Add second axis
        self.left = self.ax
        self.right = self.ax.twinx()
        self.right.axis["right"].toggle(all=True)
        
        return None
    
    def _set_colors(self, left_plot, right_plot) -> None:
        
        self.left.axis["left"].label.set_color(left_plot.get_color())
        self.right.axis["right"].label.set_color(right_plot.get_color())
        
        return None
    

class ParasiteAxis(DoubleAxis):
    
    def __init__(self, config: PlotConfig) -> None:
        
        super().__init__(config)
        
        # Add parasite axis
        self.parax = self.left.twinx()
        self.parax.axis["right"] = self.parax.new_fixed_axis(loc="right",
                                                             offset=(70, 0))
        self.parax.axis["right"].toggle(all=True)
        
    def post(self):
        
        if self.config.legend:
            self.ax.legend()
            
        # self.fig.tight_layout()
        self.fig.subplots_adjust(left=0.05, right=0.85)

        if self.config.save:
            Path(self.config.save_dir).mkdir(parents=True, exist_ok=True)
            self.path = f"{self.config.save_dir}/{self.config.save_name}"
            self.fig.savefig(self.path)

        if self.config.show:
            plt.show(block=True)

        plt.close(self.fig) 


def plot_orbit(time: np.ndarray, r: np.ndarray, u: np.ndarray) -> None:
    """Plot components and magnitude of positon and velocity vectors.

    Generates two files: components.png, which contains a plot with the
    evolution of the components of the position and velocity vectors, and
    magnitudes.png, which contains a plot with the evolution of their
    magnitudes. Files are saved to a dedicated 'plots' directory, which
    is created if it does not exist.

    :param time: Time array.
    :param r: Position array.
    :param u: Velocity array.
    """

    # Create directory for plots if it does not exist
    Path("plots").mkdir(parents=True, exist_ok=True)

    # Position and velocity components
    fig, axs = plt.subplots(figsize=(10, 6), nrows=2, ncols=1, sharex=True)
    fig.suptitle("Position and velocity components for FORMOSAT 7-6")
    axs[0].plot(time, r[:, 0], label=r"$x$")
    axs[0].plot(time, r[:, 1], label=r"$y$")
    axs[0].plot(time, r[:, 2], label=r"$z$")
    axs[0].set_ylabel("Position [m]")
    axs[0].legend()
    axs[1].plot(time, u[:, 0], label=r"$u$")
    axs[1].plot(time, u[:, 1], label=r"$v$")
    axs[1].plot(time, u[:, 2], label=r"$w$")
    axs[1].set_ylabel("Velocity [m/s]")
    axs[1].legend()
    axs[1].set_xlabel("Julian date")
    fig.savefig("plots/components.png")

    # Magnitude of position and velocity vectors
    r_mag = np.linalg.norm(r, axis=1)
    u_mag = np.linalg.norm(u, axis=1)
    fig2, axs2 = plt.subplots(figsize=(10, 6), nrows=2, ncols=1, sharex=True)
    fig2.suptitle("Magnitude of position and velocity vectors")
    axs2[0].plot(time, r_mag)
    axs2[0].set_ylabel("Radius [m]")
    axs2[1].plot(time, u_mag)
    axs2[1].set_ylabel("Velocity [m/s]")
    axs2[1].set_xlabel("Julian date")
    fig2.savefig("plots/magnitudes.png")

    return None


def plot_elements(t: Vector,
                  elements: KeplerianState, fig_name: str,
                  plot_dir: str = "plots", show: bool = False) -> str:
    """Plot evolution of orbital elements.

    :param t: Epoch in which the orbital elements are known.
    :param elements: Orbital elements.
    :param fig_name: Name of file to save the figure.
    :param plot_dir: Directory to save the figure.
    :param show: If True, show the plot.
    :returns: Path to the saved figure.
    """

    # Create directory for plots if it does not exist
    Path(plot_dir).mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(10, 6))
    fig.suptitle("Evolution of orbital elements for FORMOSAT 7-6")

    # Plot semi-major axis
    ax[0, 0].plot(t, elements.a * 1e-3)
    ax[0, 0].set_ylabel(r"$a$ [km]")
    ax[0, 0].set_xlabel("Julian day")

    # Plot eccentricity
    ax[0, 1].plot(t, elements.e)
    ax[0, 1].set_ylabel(r"$e$")
    ax[0, 1].set_xlabel("Julian day")
    ax[0, 1].set_yticks(np.arange(0.0008, 0.0015, 0.0002))

    # Plot inclination
    ax[1, 0].plot(t, elements.i)
    ax[1, 0].set_ylabel(r"$i$ [deg]")
    ax[1, 0].set_xlabel("Julian day")

    # Plot argument of periapsis
    ax[1, 1].plot(t, elements.omega)
    ax[1, 1].set_ylabel(r"$\omega$ [deg]")
    ax[1, 1].set_xlabel("Julian day")
    ax[1, 1].set_yticks(np.arange(0., 361., 90.))

    # Plot right ascension of ascending node
    ax[2, 0].plot(t, elements.Omega)
    ax[2, 0].set_ylabel(r"$\Omega$ [deg]")
    ax[2, 0].set_xlabel("Julian day")
    ax[2, 0].set_yticks(np.arange(298, 302., 1.))

    # Plot true anomaly
    ax[2, 1].plot(t, elements.nu)
    ax[2, 1].set_ylabel(r"$\theta$ [deg]")
    ax[2, 1].set_xlabel("Julian day")
    ax[2, 1].set_yticks(np.arange(0., 361., 90.))

    fig.tight_layout()
    fig.savefig(f"{plot_dir}/{fig_name}")

    if show:
        plt.show

    return f"{plot_dir}/{fig_name}"


def compare_cartesian_orbits(t: Vector, s1: CartesianState,
                             s2: CartesianState, fig_name: str,
                             plot_dir: str = "plots", log: bool = False,
                             show: bool = False) -> str:
    """Compare two cartesian orbits.

    :param t: Epochs.
    :param s1: First orbit.
    :param s2: Second orbit.
    :param fig_name: Name of the figure.
    :param plot_dir: Directory to save the figure.
    :param log: If True, plot in log scale.
    :param show: If True, show the plot.
    :returns: Path to the saved figure.
    """

    # Create directory for plots if it does not exist
    Path(plot_dir).mkdir(parents=True, exist_ok=True)
    fig_path = f"{plot_dir}/{fig_name}"

    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(10, 6))
    fig.suptitle("Normalized residuals after conversion to keplerian elements")

    delta_s = s1 - s2
    dx = delta_s.x / np.mean(s1.x)
    dy = delta_s.y / np.mean(s1.y)
    dz = delta_s.z / np.mean(s1.z)
    ddx = delta_s.dx / np.mean(s1.dx)
    ddy = delta_s.dy / np.mean(s1.dy)
    ddz = delta_s.dz / np.mean(s1.dz)

    if log:
        # Take absolute value of residuals
        dx = np.abs(dx)
        dy = np.abs(dy)
        dz = np.abs(dz)
        ddx = np.abs(ddx)
        ddy = np.abs(ddy)
        ddz = np.abs(ddz)

        # Add absolute value to labels
        axs[0, 0].set_ylabel(r"$|\Delta x\ /\ \mu_x|$")
        axs[0, 1].set_ylabel(r"$|\Delta y\ /\ \mu_y|$")
        axs[1, 0].set_ylabel(r"$|\Delta z\ /\ \mu_z|$")
        axs[1, 1].set_ylabel(r"$|\Delta \dot{x}\ /\ \mu_{\dot{x}}|$")
        axs[2, 0].set_ylabel(r"$|\Delta \dot{y}\ /\ \mu_{\dot{y}}|$")
        axs[2, 1].set_ylabel(r"$|\Delta \dot{z}\ /\ \mu_{\dot{z}}|$")
    else:
        # Define regular labels
        axs[0, 0].set_ylabel(r"$\Delta x\ /\ \mu_x$")
        axs[0, 1].set_ylabel(r"$\Delta y\ /\ \mu_y$")
        axs[1, 0].set_ylabel(r"$\Delta z\ /\ \mu_z$")
        axs[1, 1].set_ylabel(r"$\Delta \dot{x}\ /\ \mu_{\dot{x}}$")
        axs[2, 0].set_ylabel(r"$\Delta \dot{y}\ /\ \mu_{\dot{y}}$")
        axs[2, 1].set_ylabel(r"$\Delta \dot{z}\ /\ \mu_{\dot{z}}$")

    # Plot x
    axs[0, 0].plot(t, dx)
    axs[0, 0].set_xlabel("Julian day")

    # Plot y
    axs[0, 1].plot(t, dy)
    axs[0, 1].set_xlabel("Julian day")

    # Plot z
    axs[1, 0].plot(t, dz)
    axs[1, 0].set_xlabel("Julian day")

    # Plot dx
    axs[1, 1].plot(t, ddx)
    axs[1, 1].set_xlabel("Julian day")

    # Plot dy
    axs[2, 0].plot(t, ddy)
    axs[2, 0].set_xlabel("Julian day")

    # Plot dz
    axs[2, 1].plot(t, ddz)
    axs[2, 1].set_xlabel("Julian day")

    # Set yscale to log if requested
    if log:
        for ax in axs.flatten():
            ax.set_yscale("log")

    fig.tight_layout()
    fig.savefig(fig_path)

    if show:
        plt.show()

    return fig_path


def plot_spherical_coordinates(t: Vector, s: SphericalGeocentric,
                               fig_name: str, plot_dir: str = "plots",
                               show: bool = False) -> str:
    """Plot traces and orbital radius

    :param t: Time vector
    :param s: Time series of position vectors in spherical cooridinates
    :param fig_name: Name of the figure.
    :param plot_dir: Directory to save the figure.
    :param show: If True, show the plot.
    :returns: Path to the saved figure.
    """

    # Create directory for plots if it does not exist
    Path(plot_dir).mkdir(parents=True, exist_ok=True)
    fig_path = f"{plot_dir}/{fig_name}"

    # Remove discontinuities in longitude
    assert isinstance(s.lon, np.ndarray)
    diff = (s.lon[:-1] - s.lon[1:]) > 200.
    diff = np.concatenate((diff, [False]))
    s.lon[diff] = np.nan

    # Plot traces and radius
    fig, (trace_ax, radius_ax) = plt.subplots(2, 1, figsize=(10, 6))
    fig.suptitle("Orbital traces and radius for FORMOSAT 7-6")

    trace_ax.plot(s.lon, s.lat)
    trace_ax.set_xlabel("Longitude [deg]")
    trace_ax.set_ylabel("Latitude [deg]")
    trace_ax.set_xticks(np.arange(-180, 181, 60))

    radius_ax.plot(t, s.r)
    radius_ax.set_xlabel("Julian day")
    radius_ax.set_ylabel("Radius [m]")

    # Save figure and show if requested
    fig.tight_layout()
    fig.savefig(fig_path)

    if show:
        plt.show()

    return fig_path


@overload
def compare_spherical(t: Vector, s: CartesianState,
                      ref: CartesianState, config: PlotConfig) -> str: ...


@overload
def compare_spherical(t: Vector, s: SphericalGeocentric,
                      ref: SphericalGeocentric, config: PlotConfig) -> str: ...


def compare_spherical(t, s, ref, config) -> str:

    # Convert to spherical coordinates if needed and calculate residuals
    if isinstance(ref, CartesianState):
        ref = cartesian2spherical(t, ref)
        s = cartesian2spherical(t, s)
    elif isinstance(ref, SphericalGeocentric):
        pass
    else:
        raise TypeError("Input must be cartesian or spherical state vector")

    dr, dlat, dlon = s.residual(ref)

    # Plot data
    fig = plt.figure(figsize=config.figsize)
    rax = host_subplot(111,
                       axes_class=axisartist.Axes, figure=fig)  # type: ignore
    plt.subplots_adjust(right=0.8)
    lax = rax.twinx()
    lox = rax.twinx()
    lox.axis["right"] = lox.new_fixed_axis(loc="right", offset=(70, 0))

    lax.axis["right"].toggle(all=True)
    lox.axis["right"].toggle(all=True)

    if config.title is not None:
        fig.suptitle(config.title)

    if config.yscale is not None:
        lax.set_yscale(config.yscale)
        lox.set_yscale(config.yscale)
        rax.set_yscale(config.yscale)

    rad, = rax.plot(t, dr, label=r"$\Delta r$")
    lat, = lax.plot(t, dlat, label=r"$\Delta \varphi$")
    lon, = lox.plot(t, dlon, label=r"$\Delta \lambda$")

    rax.set(xlabel="Julian day", ylabel="Radial residual [m]")
    lax.set(ylabel="Latitudinal residual [deg]")
    lox.set(ylabel="Longitudinal residual [deg]")

    rax.axis["left"].label.set_color(rad.get_color())
    lax.axis["right"].label.set_color(lat.get_color())
    lox.axis["right"].label.set_color(lon.get_color())

    return __save_and_show(fig, config)


def plot_acceleration(t: Vector, s: CartesianState,
                      forces: ForceModel, config: PlotConfig) -> str:

    # Calculate acceleration along orbit
    ds = get_acceleration(t, s, forces)

    tau = (t - t[0]) * 24 * 60 / 95.24

    # Plot results
    fig, ax = plt.subplots(figsize=config.figsize)

    if config.title is not None:
        fig.suptitle(config.title)

    ax.plot(tau, ds.a_mag)
    ax.set_xlabel("Julian day")
    ax.set_ylabel(r"Acceleration [$m/s^2$]")

    if config.yscale is not None:
        ax.set_yscale(config.yscale)

    return __save_and_show(fig, config)


def plot_gdiff(t: Vector, s: dict[str, CartesianState],
               ref: CartesianState, config: PlotConfig) -> str:

    # Calculate global difference for each element of the sequence
    dr = {}
    du = {}
    for key, si in s.items():
        dri, dui = global_diff(si, ref)
        dr[key] = dri
        du[key] = dui

    # Create figure
    fig, (rax, uax) = plt.subplots(2, 1, figsize=(10, 6))
    if config.title is not None:
        fig.suptitle(config.title)

    for key in s:
        rax.plot(t, dr[key], label=key)
        uax.plot(t, du[key], label=key)

    rax.set_xlabel("Julian day")
    rax.set_ylabel(r"$\gamma^{pos} (t)$ [m]")
    rax.legend()

    uax.set_xlabel("Julian day")
    uax.set_ylabel(r"$\gamma^{vel} (t)$ [m/s]")
    uax.legend()

    if config.yscale is not None:
        rax.set_yscale(config.yscale)
        uax.set_yscale(config.yscale)

    fig.tight_layout()

    # Save and show if requested
    path = ""
    if config.save:
        Path(config.save_dir).mkdir(parents=True, exist_ok=True)
        path = f"{config.save_dir}/{config.save_name}"
        fig.savefig(path)

    if config.show:
        plt.show()

    return path


def plot_global_diff(
    t: Vector | Sequence[Vector],
    diff: tuple[Vector, Vector] | tuple[Sequence[Vector], Sequence[Vector]],
    config: PlotConfig
) -> str:
    """DEPRECATED: USE plot_gdiff instead"""

    # Cast everything to sequences
    dr, du = diff
    if (isinstance(t, np.ndarray) and isinstance(dr, np.ndarray)
            and isinstance(du, np.ndarray)):
        t = [t]
        dr = [dr]
        du = [du]

    # Create figure
    fig, (rax, uax) = plt.subplots(2, 1, figsize=(10, 6))
    if config.title is not None:
        fig.suptitle(config.title)

    for ti, dri, dui in zip(t, dr, du):

        label = f"{((ti[1] - ti[0]) * 86400.):.2f} s"
        rax.plot(ti, dri, label=label)
        uax.plot(ti, dui, label=label)

    rax.set_ylabel(r"$\gamma^{pos} (t)$ [m]")
    rax.set_xlabel("Julian day")
    rax.legend()

    uax.set_ylabel(r"$\gamma^{vel} (t)$ [m/s]")
    uax.set_xlabel("Julian day")
    uax.legend()

    if config.yscale is not None:
        rax.set_yscale(config.yscale)
        uax.set_yscale(config.yscale)

    fig.tight_layout()

    # Save figure and show if requested
    fig_path = ""
    if config.save:
        Path(config.save_dir).mkdir(parents=True, exist_ok=True)
        fig_path = f"{config.save_dir}/{config.save_name}"
        fig.savefig(fig_path)
    if config.show:
        plt.show()

    return fig_path


def plot_cglobal_diff(dt: Sequence, diff: tuple[Sequence, Sequence],
                      config: PlotConfig) -> str:

    dr, du = diff

    # Create figure
    fig, (rax, uax) = plt.subplots(2, 1, figsize=(10, 6))
    if config.title is not None:
        fig.suptitle(config.title)

    # Plot difference in position
    rax.plot(dt, dr, "o-")
    rax.set_ylabel(r"$\gamma^{(RMS, pos)} (t)$ [m]")
    rax.set_xlabel("Time step [s]")

    # Plot difference in velocity
    uax.plot(dt, du, "o-")
    uax.set_ylabel(r"$\gamma^{(RMS, vel)} (t)$ [m/s]")
    uax.set_xlabel("Time step [s]")

    if config.yscale is not None:
        rax.set_yscale(config.yscale)
        uax.set_yscale(config.yscale)

    fig.tight_layout()

    # Save figure and show if requested
    fig_path = ""
    if config.save:
        Path(config.save_dir).mkdir(parents=True, exist_ok=True)
        fig_path = f"{config.save_dir}/{config.save_name}"
        fig.savefig(fig_path)
    if config.show:
        plt.show()

    return fig_path


def plot_gdiff_position(t: Sequence[Vector], dr: Sequence[Vector],
                        config: PlotConfig) -> str:

    # Check input to ensure we get 6 components
    if (len(t) != 6) or (len(dr) != 6):
        raise ValueError("Unexpected length of time or dr sequence")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    if config.title is not None:
        fig.suptitle(config.title)

    # Plot difference in position
    labels = [r"$\delta x$", r"$\delta y$", r"$\delta z$",
              r"$\delta \dot x$", r"$\delta \dot y$", r"$\delta \dot z$"]

    for ti, dri, label in zip(t, dr, labels):
        ax.plot(ti, dri * 1e-3, label=label)

    ax.set_xlabel("Julian date")
    ax.set_ylabel(r"$\gamma^{pos}(t) [km]$")
    ax.legend()

    if config.yscale is not None:
        ax.set_yscale(config.yscale)
    fig.tight_layout()

    # Save figure and show if requested
    fig_path = ""
    if config.save:
        Path(config.save_dir).mkdir(parents=True, exist_ok=True)
        fig_path = f"{config.save_dir}/{config.save_name}"
        fig.savefig(fig_path)
    if config.show:
        plt.show()

    return fig_path


def plot_cgdiff_euler_ds(t, cgdiff, config):

    fig, (rax, uax) = plt.subplots(2, 1, figsize=(10, 6))

    if config.title is not None:
        fig.suptitle(config.title)

    # Express perturbation as percentage
    t = np.array(t) * 100.
    dx, dy, dz, ddx, ddy, ddz = cgdiff

    # Plot position differences
    rax.plot(t, dx, "o-", label=r"$\delta x$")
    rax.plot(t, dy, "o-", label=r"$\delta y$")
    rax.plot(t, dz, "o-", label=r"$\delta z$")
    rax.set_xlabel("Percentage of perturbation")
    rax.set_ylabel(r"$\gamma^{(RMS, pos)}(t)$ [m]")
    rax.legend()

    # Plot velocity differences
    uax.plot(t, ddx, "o-", label=r"$\delta \dot x$")
    uax.plot(t, ddy, "o-", label=r"$\delta \dot y$")
    uax.plot(t, ddz, "o-", label=r"$\delta \dot z$")
    uax.set_xlabel("Percentage of perturbation")
    uax.set_ylabel(r"$\gamma^{(RMS, pos)}(t)$ [m]")
    uax.legend()

    fig.tight_layout()

    if config.yscale is not None:
        rax.set_yscale(config.yscale)
        uax.set_yscale(config.yscale)
    rax.set_xscale("log")
    uax.set_xscale("log")

    fig_path = ""
    if config.save:
        Path(config.save_dir).mkdir(parents=True, exist_ok=True)
        fig_path = f"{config.save_dir}/{config.save_name}"
        fig.savefig(fig_path)
    if config.show:
        plt.show()

    return fig_path
