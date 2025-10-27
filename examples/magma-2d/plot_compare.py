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
    """Command line interface.
    """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--no-gui", action="store_false", dest="show_plot", default=True, help="Do not display plot.")

    args = parser.parse_args()
    PlotApp().run(
            show_plot=args.show_plot,
        )

class PlotApp:
    """Application to plot displacement and fluid pressure."""
    OUTPUT_DIR = pathlib.Path("output")
    FILENAME_ORIG = OUTPUT_DIR / "step01_inflation-domain.h5"
    FILENAME_ADAPT = OUTPUT_DIR / "step01b_inflation-domain.h5"

    N_ROWS = 2
    N_COLS = 1
    STYLE_UNIFORM = "k-o"
    STYLE_ADAPT = "r--x"

    def run(self, show_plot: bool):
        orig = self._load_data(self.FILENAME_ORIG)
        adapt = self._load_data(self.FILENAME_ADAPT)

        self._plot(orig=orig, adapt=adapt, show_plot=show_plot)

    def _load_data(self, filename: pathlib.Path):
        """Load simulation results."""
        h5 = h5py.File(filename)
        data = {
            "time": h5["/time"][:].squeeze(),
            "vertices": h5["/geometry/vertices"][:],
            "displacement": h5["/vertex_fields/displacement"][:],
            "fluid_pressure": h5["/vertex_fields/pressure"][:]
        }
        h5.close()
        return data

    def _plot(self, orig: dict, adapt: dict, show_plot: bool):
        figure = pyplot.figure(figsize=(7.0, 6.0), layout="tight")

        # Plot displacement
        axes = self._setup_axes(figure, index=1)
        i_loc = numpy.argmax(numpy.abs(orig["displacement"][-1,:,1]))
        axes.plot(orig["time"]/year.value, orig["displacement"][:,i_loc,1], "k-o")
        axes.plot(adapt["time"]/year.value, adapt["displacement"][:,i_loc,1], "r--x")
        vertex = orig["vertices"][i_loc]
        axes.set_title(f"Vertical displacement at x={vertex[0]:.0f}, y={vertex[1]:.0f}")
        axes.set_ylabel("Y displacement, m")

        # Plot fluid pressure
        axes = self._setup_axes(figure, index=2)
        orig_p = orig["fluid_pressure"][:,:,0]
        i_loc = numpy.argmax(numpy.abs(orig_p[-1,:] - orig_p[1,:]))
        axes.plot(orig["time"]/year.value, orig["fluid_pressure"][:,i_loc,0]/MPa.value, "k-o")
        axes.plot(adapt["time"]/year.value, adapt["fluid_pressure"][:,i_loc,0]/MPa.value, "r--x")
        vertex = orig["vertices"][i_loc]
        axes.set_title(f"Fluid pressure at x={vertex[0]:.0f}, y={vertex[1]:.0f}")
        axes.set_ylabel("Fluid pressure, MPa")

        if show_plot:
            pyplot.show()
        figure.savefig("step01-compare.pdf")

    def _setup_axes(self, figure, index: int):
        axes = figure.add_subplot(self.N_ROWS, self.N_COLS, index)
        axes.spines[["top", "right"]].set_visible(False)
        axes.set_xlabel("Time, year")
        return axes

if __name__ == "__main__":
    cli()