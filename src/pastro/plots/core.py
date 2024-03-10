from dataclasses import dataclass
import matplotlib.pyplot as plt
from pathlib import Path
from ..types import Vector
import numpy as np
from typing import TypeVar


# TODO: How do I add labels if more than one line in a plot
# TODO: Optional argument in plot method to prevent postprocessing
@dataclass
class PlotSetup:
    """Plot configuration.

    :param figsize: Size of the figure [(10, 6)]
    :param subplots: Number of subplots as rows and columns [(1, 1)]
    :param right_axis: If True, show right axis [False]
    :param parasite_axis: If True, show parasite axis [False]
    :param top_axis: If True, show top axis [False]
    :param title: Title of the figure [None]
    :param xlabel: Label for the x-axis [None]
    :param tlabel: Label for the top x-axis [None]
    :param ylabel: Label for the y-axis [None]
    :param rlabel: Label for the right y-axis [None]
    :param plabel: Label for the parasite y-axis [None]
    :param legend: If True, show legend [False]
    :param xscale: Scale for the x-axis [None]
    :param tscale: Scale for the top x-axis [None]
    :param yscale: Scale for the y-axis [None]
    :param rscale: Scale for the right y-axis [None]
    :param pscale: Scale for the parasite y-axis [None]
    :param xlim: Limits for the x-axis [None]
    :param tlim: Limits for the top x-axis [None]
    :param ylim: Limits for the y-axis [None]
    :param rlim: Limits for the right y-axis [None]
    :param plim: Limits for the parasite y-axis [None]
    :param show: If True, show the plot [True]
    :param save: If True, save the plot [False]
    :param dir: Directory to save the plot ["plots"]
    :param name: Name of the file to save with extension [""]
    """

    # Basic figure configuration
    figsize: tuple[float, float] = (10, 6)
    subplots: tuple[int, int] = (1, 1)

    # Selection of alternative axes
    right_axis: bool = False
    parasite_axis: bool = False
    top_axis: bool = False

    # Title and labels
    title: str | None = None
    xlabel: str | None = None
    tlabel: str | None = None
    ylabel: str | None = None
    rlabel: str | None = None
    plabel: str | None = None
    legend: bool = False

    # Scales
    xscale: str | None = None
    tscale: str | None = None
    yscale: str | None = None
    rscale: str | None = None
    pscale: str | None = None

    # Limits
    xlim: tuple[float, float] | None = None
    tlim: tuple[float, float] | None = None
    ylim: tuple[float, float] | None = None
    rlim: tuple[float, float] | None = None
    plim: tuple[float, float] | None = None

    # Save and show configuration
    show: bool = True
    save: bool = False
    dir: str = "plots"
    name: str = ""


T = TypeVar("T", bound="BasePlot")


class BasePlot:

    def __init__(self, setup: PlotSetup = PlotSetup(), _fig=None, _ax=None) -> None:

        self.COLOR_CYCLE = iter(
            [
                "#1f77b4",
                "#aec7e8",
                "#ff7f0e",
                "#ffbb78",
                "#2ca02c",
                "#98df8a",
                "#d62728",
                "#ff9896",
                "#9467bd",
                "#c5b0d5",
                "#8c564b",
                "#c49c94",
                "#e377c2",
                "#f7b6d2",
                "#7f7f7f",
                "#c7c7c7",
                "#bcbd22",
                "#dbdb8d",
                "#17becf",
                "#9edae5",
            ]
        )

        self.setup: PlotSetup = setup

        # Create figure and left axis
        if _fig is None and _ax is None:
            self.fig, self.ax = plt.subplots(figsize=self.setup.figsize)
        elif _fig is not None and _ax is not None:
            self.fig = _fig
            self.ax = _ax
        else:
            raise ValueError("Provide both figure and axis or none of them")

        # Set title and labels
        if self.setup.title is not None:
            self.fig.suptitle(self.setup.title)

        if self.setup.xlabel is not None:
            self.ax.set_xlabel(self.setup.xlabel)
        if self.setup.xscale is not None:
            self.ax.set_xscale(self.setup.xscale)
        if self.setup.xlim is not None:
            self.ax.set_xlim(self.setup.xlim)

        if self.setup.ylabel is not None:
            self.ax.set_ylabel(self.setup.ylabel)
        if self.setup.yscale is not None:
            self.ax.set_yscale(self.setup.yscale)
        if self.setup.ylim is not None:
            self.ax.set_ylim(self.setup.ylim)

        self.path: str = ""

        return None

    def _formatter(self, x, pos):

        if x == 0.0:
            return f"{x:.0f}"
        elif (np.abs(x) > 0.01) and (np.abs(x) <= 1e3):
            return f"{x}"
        else:
            a, b = f"{x:.0e}".split("e")
            bsign = "-" if a[0] == "-" else ""
            esign = "-" if b[0] == "-" else ""
            exp = int(b[1:])
            n = int(a[0]) if bsign == "" else int(a[1])
            return f"${bsign}{n}" + r"\cdot 10^{" + f"{esign}{exp}" + r"}$"

    def _postprocess(self) -> None:

        for line in self.ax.get_lines():
            line.set_color(next(self.COLOR_CYCLE))

        return None

    def postprocess(self) -> None:

        self._postprocess()

        if self.setup.legend:
            self.ax.legend()

        if self.setup.save:
            Path(self.setup.dir).mkdir(parents=True, exist_ok=True)
            self.path = f"{self.setup.dir}/{self.setup.name}"
            self.fig.savefig(self.path)

        if self.setup.show:
            plt.show(block=True)
            plt.close(self.fig)

        return None

    def _plot(
        self, x: Vector, y: Vector, label: str | None = None, axis: str | None = None
    ) -> None:
        raise NotImplementedError

    def add_line(
        self, x: Vector, y: Vector, label: str | None = None, axis: str | None = None
    ) -> None:
        self._plot(x, y, label=label, axis=axis)
        return None

    def __call__(self) -> str:
        self.postprocess()
        return self.path


class SinglePlot(BasePlot):

    def _plot(
        self, x: Vector, y: Vector, label: str | None = None, axis: str | None = "left"
    ) -> None:

        self.ax.plot(x, y, label=label)
        return None

    def plot(self, x: Vector, y: Vector) -> str:

        self._plot(x, y)
        return self.__call__()


# class SinglePlot:

#     def __init__(
#         self, setup: PlotSetup = PlotSetup(), _figure=None, _axis=None
#     ) -> None:

#         self.COLOR_CYCLE = iter(
#             [
#                 "#1f77b4",
#                 "#aec7e8",
#                 "#ff7f0e",
#                 "#ffbb78",
#                 "#2ca02c",
#                 "#98df8a",
#                 "#d62728",
#                 "#ff9896",
#                 "#9467bd",
#                 "#c5b0d5",
#                 "#8c564b",
#                 "#c49c94",
#                 "#e377c2",
#                 "#f7b6d2",
#                 "#7f7f7f",
#                 "#c7c7c7",
#                 "#bcbd22",
#                 "#dbdb8d",
#                 "#17becf",
#                 "#9edae5",
#             ]
#         )

#         self.setup: PlotSetup = setup

#         # Create figure and left axis
#         if _figure is None and _axis is None:
#             self.fig, self.ax = plt.subplots(figsize=self.setup.figsize)
#         elif _figure is not None and _axis is not None:
#             self.fig = _figure
#             self.ax = _axis
#         else:
#             raise ValueError("Provide both figure and axis or none of them")

#         # Set title and labels
#         if self.setup.title is not None:
#             self.fig.suptitle(self.setup.title)

#         if self.setup.xlabel is not None:
#             self.ax.set_xlabel(self.setup.xlabel)
#         if self.setup.xscale is not None:
#             self.ax.set_xscale(self.setup.xscale)
#         if self.setup.xlim is not None:
#             self.ax.set_xlim(self.setup.xlim)

#         if self.setup.ylabel is not None:
#             self.ax.set_ylabel(self.setup.ylabel)
#         if self.setup.yscale is not None:
#             self.ax.set_yscale(self.setup.yscale)
#         if self.setup.ylim is not None:
#             self.ax.set_ylim(self.setup.ylim)

#         self.path: str = ""

#         return None

#     def _formatter(self, x, pos):

#         if x == 0.0:
#             return f"{x:.0f}"
#         elif (np.abs(x) > 0.01) and (np.abs(x) < 1e4):
#             return f"{x}"
#         else:
#             a, b = f"{x:.0e}".split("e")
#             bsign = "-" if a[0] == "-" else ""
#             esign = "-" if b[0] == "-" else ""
#             exp = int(b[1:])
#             n = int(a[0]) if bsign == "" else int(a[1])
#             return f"${bsign}{n}" + r"\cdot 10^{" + f"{esign}{exp}" + r"}$"

#     def _postprocess(self) -> None:

#         for line in self.ax.get_lines():
#             line.set_color(next(self.COLOR_CYCLE))

#         return None

#     def postprocess(self) -> None:

#         self._postprocess()

#         if self.setup.legend:
#             self.ax.legend()

#         if self.setup.save:
#             Path(self.setup.dir).mkdir(parents=True, exist_ok=True)
#             self.path = f"{self.setup.dir}/{self.setup.name}"
#             self.fig.savefig(self.path)

#         if self.setup.show:
#             plt.show(block=True)
#             plt.close(self.fig)

#         return None

#     def add_line(self, x: Vector, y: Vector, label: str | None = None) -> None:

#         self._plot(x, y, label=label)
#         return None

#     def plot(self, x: Vector, y: Vector, label: str | None = None) -> str | None:

#         self._plot(x, y, label=label)
#         self.postprocess()
#         return self.path

#     def _plot(self, x: Vector, y: Vector, label: str | None = None) -> None:

#         self.ax.plot(x, y, label=label)
#         return None


class DoubleAxis(BasePlot):

    def __init__(self, setup: PlotSetup = PlotSetup(), _fig=None, _ax=None) -> None:

        setup.right_axis = True

        super().__init__(setup, _fig, _ax)

        self.left = self.ax
        self.right = self.left.twinx()

        if self.setup.rlabel is not None:
            self.right.set_ylabel(self.setup.rlabel)
        if self.setup.rscale is not None:
            self.right.set_yscale(self.setup.rscale)
        if self.setup.rlim is not None:
            self.right.set_ylim(self.setup.rlim)

        return None

    def _postprocess(self) -> None:

        super()._postprocess()

        for line in self.right.get_lines():
            line.set_color(next(self.COLOR_CYCLE))
        right_line = self.right.get_lines()[0]
        self.right.yaxis.label.set_color(right_line.get_color())
        self.fig.subplots_adjust(left=0.1, right=0.8)
        self.right.yaxis.set_major_formatter(self._formatter)

        return None

    def _plot(
        self, x: Vector, y: Vector, label: str | None = None, axis: str | None = None
    ) -> None:

        if axis == "left":
            self.ax.plot(x, y, label=label)
        elif axis == "right":
            self.right.plot(x, y, label=label)
        else:
            raise ValueError("Axis must be either 'left' or 'right'")
        return None

    def plot(self, x: Vector, y_left: Vector, y_right: Vector) -> str:

        self._plot(x, y_left, axis="left")
        self._plot(x, y_right, axis="right")
        return self.__call__()


# class DoubleAxis(SinglePlot):

#     def __init__(
#         self, setup: PlotSetup = PlotSetup(), _figure=None, _axis=None
#     ) -> None:

#         super().__init__(setup, _figure, _axis)

#         if not self.setup.right_axis:
#             raise ValueError("Right axis must be enabled to use DoubleAxis")

#         self.left = self.ax
#         self.right = self.left.twinx()

#         if self.setup.rlabel is not None:
#             self.right.set_ylabel(self.setup.rlabel)
#         if self.setup.rscale is not None:
#             self.right.set_yscale(self.setup.rscale)
#         if self.setup.rlim is not None:
#             self.right.set_ylim(self.setup.rlim)

#         return None

#     def _postprocess(self) -> None:

#         super()._postprocess()

#         for line in self.right.get_lines():
#             line.set_color(next(self.COLOR_CYCLE))
#         right_line = self.right.get_lines()[0]
#         self.right.yaxis.label.set_color(right_line.get_color())
#         self.fig.subplots_adjust(left=0.1, right=0.8)
#         self.right.yaxis.set_major_formatter(self._formatter)

#         return None

#     def _plot(self, axis: str, x: Vector, y: Vector, label: str | None = None) -> None:

#         if axis == "left":
#             self.ax.plot(x, y, label=label)
#         elif axis == "right":
#             self.right.plot(x, y, label=label)
#         else:
#             raise ValueError("Axis must be either 'left' or 'right'")
#         return None

#     def add_line(
#         self, axis: str, x: Vector, y: Vector, label: str | None = None
#     ) -> None:

#         self._plot(axis, x, y, label=label)
#         return None

#     def plot(self, x: Vector, y: Vector, x_right: Vector, y_right: Vector) -> str:

#         self._plot("left", x, y)
#         self._plot("right", x_right, y_right)
#         self.postprocess()
#         return self.path

#     # def __plot(self, x: Vector, y: Vector, x_right: Vector, y_right: Vector) -> None:

#     #     self.ax.plot(x, y)
#     #     self.right.plot(x_right, y_right)
#     #     return None

#     # def plot(self, *args, **kwargs) -> str | None:

#     #     self.__plot(*args, **kwargs)
#     #     self.postprocess()
#     #     return self.path


class ParasiteAxis(BasePlot):

    def __init__(self, setup: PlotSetup = PlotSetup(), _fig=None, _ax=None) -> None:

        setup.right_axis = True
        setup.parasite_axis = True

        super().__init__(setup, _fig, _ax)

        self.left = self.ax
        self.right = self.left.twinx()

        if self.setup.rlabel is not None:
            self.right.set_ylabel(self.setup.rlabel)
        if self.setup.rscale is not None:
            self.right.set_yscale(self.setup.rscale)
        if self.setup.rlim is not None:
            self.right.set_ylim(self.setup.rlim)

        self.parax = self.left.twinx()
        self.parax.spines.right.set_position(("axes", 1.13))

        if self.setup.plabel is not None:
            self.parax.set_ylabel(self.setup.plabel)
        if self.setup.pscale is not None:
            self.parax.set_yscale(self.setup.pscale)
        if self.setup.plim is not None:
            self.parax.set_ylim(self.setup.plim)

        return None

    def _postprocess(self) -> None:

        super()._postprocess()

        for line in self.right.get_lines():
            line.set_color(next(self.COLOR_CYCLE))
        right_line = self.right.get_lines()[0]
        self.right.yaxis.label.set_color(right_line.get_color())
        self.fig.subplots_adjust(left=0.1, right=0.8)
        self.right.yaxis.set_major_formatter(self._formatter)

        for line in self.parax.get_lines():
            line.set_color(next(self.COLOR_CYCLE))
        parax_line = self.parax.get_lines()[0]
        self.parax.yaxis.label.set_color(parax_line.get_color())
        self.parax.yaxis.set_major_formatter(self._formatter)

        return None

    def _plot(
        self, x: Vector, y: Vector, label: str | None = None, axis: str | None = None
    ) -> None:

        if axis == "left":
            self.ax.plot(x, y, label=label)
        elif axis == "right":
            self.right.plot(x, y, label=label)
        elif axis == "parax":
            self.parax.plot(x, y, label=label)
        else:
            raise ValueError("Axis must be either 'left', 'right' or 'parax'")
        return None

    def plot(self, x: Vector, y_left: Vector, y_right: Vector, y_parax: Vector) -> str:

        self._plot(x, y_left, axis="left")
        self._plot(x, y_right, axis="right")
        self._plot(x, y_parax, axis="parax")
        return self.__call__()


# class ParasiteAxis(DoubleAxis):

#     def __init__(
#         self, setup: PlotSetup = PlotSetup(), _figure=None, _axis=None
#     ) -> None:

#         super().__init__(setup, _figure, _axis)

#         if not self.setup.parasite_axis:
#             raise ValueError("Parasite axis must be enabled to use ParasiteAxis")

#         self.parax = self.left.twinx()
#         self.parax.spines.right.set_position(("axes", 1.13))

#         if self.setup.plabel is not None:
#             self.parax.set_ylabel(self.setup.plabel)
#         if self.setup.pscale is not None:
#             self.parax.set_yscale(self.setup.pscale)
#         if self.setup.plim is not None:
#             self.parax.set_ylim(self.setup.plim)

#         return None

#     def _postprocess(self) -> None:

#         super()._postprocess()

#         for line in self.parax.get_lines():
#             line.set_color(next(self.COLOR_CYCLE))
#         parax_line = self.parax.get_lines()[0]
#         self.parax.yaxis.label.set_color(parax_line.get_color())
#         self.parax.yaxis.set_major_formatter(self._formatter)

#         return None

#     def __plot(
#         self,
#         x: Vector,
#         y: Vector,
#         x_right: Vector,
#         y_right: Vector,
#         x_parax: Vector,
#         y_parax: Vector,
#     ) -> None:

#         self.ax.plot(x, y)
#         self.right.plot(x_right, y_right)
#         self.parax.plot(x_parax, y_parax)
#         return None

#     def plot(self, *args, **kwargs) -> str | None:

#         self.__plot(*args, **kwargs)
#         self.postprocess()
#         return self.path


class MultiPlot:

    def __init__(self, setup: PlotSetup) -> None:

        self.setup = setup

        self.rows, self.cols = self.setup.subplots
        self.fig, self.axes = plt.subplots(
            self.rows, self.cols, figsize=self.setup.figsize, layout="tight"
        )

        self.ax_list = iter(self.axes.ravel())

        if self.setup.title is not None:
            self.fig.suptitle(self.setup.title)

        self.path = ""

        return None

    def add_plot(self, setup: PlotSetup = PlotSetup(), type: type[T] = SinglePlot) -> T:
        setup.show = False
        setup.save = False
        return type(setup, _fig=self.fig, _ax=next(self.ax_list))

    def __call__(self) -> str:

        if self.setup.save:
            Path(self.setup.dir).mkdir(parents=True, exist_ok=True)
            self.path = f"{self.setup.dir}/{self.setup.name}"
            self.fig.savefig(self.path)

        if self.setup.show:
            plt.show(block=True)
            plt.close(self.fig)

        return self.path


class GenericPlot:

    def __init__(
        self, setup: PlotSetup = PlotSetup(), _figure=None, _axis=None
    ) -> None:

        self.COLOR_CYCLE = iter(
            [
                "#1f77b4",
                "#aec7e8",
                "#ff7f0e",
                "#ffbb78",
                "#2ca02c",
                "#98df8a",
                "#d62728",
                "#ff9896",
                "#9467bd",
                "#c5b0d5",
                "#8c564b",
                "#c49c94",
                "#e377c2",
                "#f7b6d2",
                "#7f7f7f",
                "#c7c7c7",
                "#bcbd22",
                "#dbdb8d",
                "#17becf",
                "#9edae5",
            ]
        )

        self.setup = setup

        # Create figure and axis
        if _figure is None and _axis is None:
            self.fig, self.ax = plt.subplots(figsize=self.setup.figsize)
        elif _figure is not None and _axis is not None:
            self.fig = _figure
            self.ax = _axis
        else:
            raise ValueError("Provide both figure and axis or none of them")

        # Set title and labels
        if self.setup.title is not None:
            self.fig.suptitle(self.setup.title)

        if self.setup.xlabel is not None:
            self.ax.set_xlabel(self.setup.xlabel)
        if self.setup.xscale is not None:
            self.ax.set_xscale(self.setup.xscale)
        if self.setup.xlim is not None:
            self.ax.set_xlim(self.setup.xlim)

        if self.setup.ylabel is not None:
            self.ax.set_ylabel(self.setup.ylabel)
        if self.setup.yscale is not None:
            self.ax.set_yscale(self.setup.yscale)
        if self.setup.ylim is not None:
            self.ax.set_ylim(self.setup.ylim)

        # Add additional axes if requested
        if self.setup.right_axis:

            self.left = self.ax
            self.right = self.left.twinx()

            if self.setup.rlabel is not None:
                self.right.set_ylabel(self.setup.rlabel)
            if self.setup.rscale is not None:
                self.right.set_yscale(self.setup.rscale)
            if self.setup.rlim is not None:
                self.right.set_ylim(self.setup.rlim)

            if self.setup.parasite_axis:
                self.parax = self.left.twinx()
                self.parax.spines.right.set_position(("axes", 1.13))

                if self.setup.plabel is not None:
                    self.parax.set_ylabel(self.setup.plabel)
                if self.setup.pscale is not None:
                    self.parax.set_yscale(self.setup.pscale)
                if self.setup.plim is not None:
                    self.parax.set_ylim(self.setup.plim)

        if self.setup.top_axis:

            self.down = self.ax
            self.up = self.down.twiny()

            if self.setup.tlabel is not None:
                self.up.set_xlabel(self.setup.tlabel)
            if self.setup.tscale is not None:
                self.up.set_xscale(self.setup.tscale)
            if self.setup.tlim is not None:
                self.up.set_xlim(self.setup.tlim)

        self.path = ""

        return None

    def _formatter(self, x, pos):

        if x == 0.0:
            return f"{x:.0f}"
        elif (np.abs(x) > 0.01) and (np.abs(x) < 1e4):
            return f"{x}"
        else:
            a, b = f"{x:.0e}".split("e")
            bsign = "-" if a[0] == "-" else ""
            esign = "-" if b[0] == "-" else ""
            exp = int(b[1:])
            n = int(a[0]) if bsign == "" else int(a[1])
            return f"${bsign}{n}" + r"\cdot 10^{" + f"{esign}{exp}" + r"}$"

    def _post(self) -> None:

        for line in self.ax.get_lines():
            line.set_color(next(self.COLOR_CYCLE))

        if self.setup.right_axis:
            for line in self.right.get_lines():
                line.set_color(next(self.COLOR_CYCLE))
            right_line = self.right.get_lines()[0]
            self.right.yaxis.label.set_color(right_line.get_color())
            self.fig.subplots_adjust(left=0.1, right=0.8)
            self.right.yaxis.set_major_formatter(self._formatter)

        if self.setup.parasite_axis:
            for line in self.parax.get_lines():
                line.set_color(next(self.COLOR_CYCLE))
            parax_line = self.parax.get_lines()[0]
            self.parax.yaxis.label.set_color(parax_line.get_color())
            self.parax.yaxis.set_major_formatter(self._formatter)

        if self.setup.top_axis:
            for line in self.up.get_lines():
                line.set_color(next(self.COLOR_CYCLE))
            up_line = self.up.get_lines()[0]
            self.up.xaxis.label.set_color(up_line.get_color())

        if self.setup.legend:
            self.ax.legend()

        if self.setup.save:
            Path(self.setup.dir).mkdir(parents=True, exist_ok=True)
            self.path = f"{self.setup.dir}/{self.setup.name}"
            self.fig.savefig(self.path)

        if self.setup.show:
            plt.show(block=True)
            plt.close(self.fig)

        return None

    def __plot(self, x: Vector, y: Vector) -> None:

        self.ax.plot(x, y)
        return None

    def plot(self, *args, **kwargs) -> str | None:

        self.__plot(*args, **kwargs)
        self._post()
        return self.path


# class MultiPlot:

#     def __init__(
#         self,
#         figure_setup: PlotSetup,
#         subplot_setups: list[PlotSetup],
#         generators: type[GenericPlot] | list[type[GenericPlot]],
#     ) -> None:

#         self.setup = figure_setup

#         self.rows, self.cols = self.setup.subplots
#         self.fig, self.axes = plt.subplots(
#             self.rows, self.cols, figsize=self.setup.figsize, layout="tight"
#         )

#         if self.setup.title is not None:
#             self.fig.suptitle(self.setup.title)

#         for setup in subplot_setups:
#             setup.show = False
#             setup.save = False

#         if isinstance(generators, type):
#             _generators = [generators for _ in range(len(subplot_setups))]
#         elif isinstance(generators, list):
#             _generators = generators
#         else:
#             raise TypeError("generators must be a type or a list of types")

#         self.generators = [
#             gen(setup, _figure=self.fig, _axis=ax)
#             for gen, setup, ax in zip(_generators, subplot_setups, self.axes.ravel())
#         ]

#         self.subplot_counter = 0

#         return None

#     def plot(self, *args, **kwargs) -> None:

#         self.generators[self.subplot_counter].plot(*args, **kwargs)
#         self.subplot_counter += 1

#         if self.subplot_counter == len(self.generators):
#             self._post()

#     def _post(self) -> None:

#         if self.setup.save:
#             Path(self.setup.dir).mkdir(parents=True, exist_ok=True)
#             self.path = f"{self.setup.dir}/{self.setup.name}"
#             self.fig.savefig(self.path)

#         if self.setup.show:
#             plt.show(block=True)
#             plt.close(self.fig)

#         return None
