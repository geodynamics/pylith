#!/usr/bin/env nemesis
"""
This script generates a plot comparing the displacement and pressure
time histories for uniform and adaptive time stepping in Step 1.
"""

# Standard Python modules
import pathlib

# Python modules included in PyLith binary installation
from matplotlib import pyplot
import h5py
import numpy

from pythia.pyre.units.time import year
from pythia.pyre.units.pressure import MPa


def cli():
    """Command line interface."""
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--step",
        action="store",
        dest="sim_step",
        default="7",
        choices=("7", "8"),
        help="Simulation step to plot.",
    )
    parser.add_argument(
        "--no-gui",
        action="store_false",
        dest="show_plot",
        default=True,
        help="Do not display plot on screen.",
    )

    args = parser.parse_args()
    PlotApp().run(
        sim_step=args.sim_step,
        show_plot=args.show_plot,
    )


class PlotApp:
    """Application to plot displacement and fluid pressure."""

    OUTPUT_DIR = pathlib.Path("output")

    N_ROWS = 2
    N_COLS = 1
    STYLE_UNIFORM = "k-o"
    STYLE_ADAPT = "r--x"
    TARGET = (-15.0e3, -25.0e3)

    def run(self, sim_step: str, show_plot: bool):
        """Run application."""
        orig = self._load_data(sim_step)
        adapt = self._load_data(sim_step + "b")

        self._plot(sim_step=sim_step, orig=orig, adapt=adapt, show_plot=show_plot)

    def _load_data(self, sim_id: str):
        """Load simulation results."""
        if sim_id[0] == "7":
            filename = self.OUTPUT_DIR / f"step0{sim_id}_twofaults_maxwell-slab.h5"
        elif sim_id[0] == "8":
            filename = self.OUTPUT_DIR / f"step0{sim_id}_twofaults_powerlaw-slab.h5"
        else:
            raise ValueError(f"Unknown simulation id '{sim_id}'.")

        h5 = h5py.File(filename)
        vertices = h5["/geometry/vertices"][:]
        data = {
            "time": h5["/time"][:].squeeze(),
            "vertices": vertices,
            "displacement": h5["/vertex_fields/displacement"][:],
            "cauchy_stress": h5["/vertex_fields/cauchy_stress"][:],
        }
        h5.close()
        return data

    def _plot(self, sim_step: str, orig: dict, adapt: dict, show_plot: bool):
        """Plot distplacement and stress time histories."""
        figure = pyplot.figure(figsize=(7.0, 6.0), layout="tight")

        # Plot displacement
        axes = self._setup_axes(figure, index=1)
        i_loc = self._find_loc(orig["vertices"], x=self.TARGET[0], y=self.TARGET[1])
        axes.plot(orig["time"] / year.value, orig["displacement"][:, i_loc, 1], "k-o")
        axes.plot(
            adapt["time"] / year.value, adapt["displacement"][:, i_loc, 1], "r--x"
        )
        vertex = orig["vertices"][i_loc]
        axes.set_title(f"Vertical displacement at x={vertex[0]:.0f}, y={vertex[1]:.0f}")
        axes.set_ylabel("Y displacement, m")

        # Plot shear stress
        axes = self._setup_axes(figure, index=2)
        i_loc = self._find_loc(orig["vertices"], x=self.TARGET[0], y=self.TARGET[1])
        axes.plot(
            orig["time"] / year.value,
            orig["cauchy_stress"][:, i_loc, 3] / MPa.value,
            "k-o",
        )
        axes.plot(
            adapt["time"] / year.value,
            adapt["cauchy_stress"][:, i_loc, 3] / MPa.value,
            "r--x",
        )
        vertex = orig["vertices"][i_loc]
        axes.set_title(f"Shear stress at x={vertex[0]:.0f}, y={vertex[1]:.0f}")
        axes.set_ylabel("Shear stress, MPa")

        if show_plot:
            pyplot.show()
        figure.savefig(f"step0{sim_step}-compare.pdf")

    def _setup_axes(self, figure, index: int):
        """Setup of axes for plotting."""
        axes = figure.add_subplot(self.N_ROWS, self.N_COLS, index)
        axes.spines[["top", "right"]].set_visible(False)
        axes.set_xlabel("Time, year")
        return axes

    def _find_loc(self, points, x, y):
        """Find index of target location given by x and y."""
        dist = numpy.sqrt(
            numpy.square(points[:, 0] - x) + numpy.square(points[:, 1] - y)
        )
        return numpy.argmin(dist)


if __name__ == "__main__":
    cli()
