#!/usr/bin/env nemesis
"""
This script generates a plot comparing the predicted solutions to
the true solution for Step 6.
"""

# Standard Python modules
import pathlib
import re

# Python modules included in PyLith binary installation
import numpy
import h5py
import matplotlib.pyplot as pyplot
from pythia.pyre.units.length import km

OUTPUT_DIR = pathlib.Path("output")

FILENAME_MODELS = OUTPUT_DIR / "step06_inversion-results.txt"
FILENAME_PRESCRIBED = OUTPUT_DIR / "step04_varslip-fault.h5"
FILENAME_RAW = OUTPUT_DIR / "slip_variable.txt"


def cli():
    """Command line interface.
    """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--models", action="store", dest="filename_models", type=str, default=FILENAME_MODELS, help="Filename of output from inversion.")
    parser.add_argument("--prescribed-slip", action="store", dest="filename_prescribed", type=str, default=FILENAME_PRESCRIBED, help="Name of HDF5 file for fault from prescribed slip simulation.")
    parser.add_argument("--raw", action="store", dest="filename_raw", type=str, default=FILENAME_RAW, help="Name of ASCII file with raw slip profile.")
    parser.add_argument("--prescribed-only", action="store_true", dest="prescribed_only", default=False, help="Show only observed slip.")
    parser.add_argument("--no-gui", action="store_false", dest="show_plot", default=True, help="Do not display plot.")

    args = parser.parse_args()
    PlotApp().run(
        filename_models=args.filename_models,
        filename_prescribed=args.filename_prescribed,
        filename_raw=args.filename_raw,
        prescribed_only=args.prescribed_only,
        show_plot=args.show_plot,
    )


class PlotApp:
    """Application to plot simple inversion results.
    """
    COLOR_RAW = "gray"
    COLOR_OBSERVED = "black"
    COLORS_MODEL = ["orange", "purple", "red", "cyan"]
    XLIM = (-25.0, 35.0)

    def run(self, filename_models: str, filename_prescribed: str, filename_raw: str, prescribed_only: bool, show_plot: bool):
        """Run plotting application.

        Args:
            filename_models: Filename of output from inversion.
            filename_prescribed: Name of HDF5 file with station output with fake observations.
            filename_raw: Name of ASCII file with raw slip values.
            prescribed_only: Show only prescribed slip.
            show_plot: Show plot window.
        """
        self._load_observed(filename_prescribed)
        self._load_raw(filename_raw)
        self._load_models(filename_models)
        self._plot(prescribed_only, show_plot)

    def _load_observed(self, filename: str):
        h5 = h5py.File(filename, "r")
        observed_coords = h5['geometry/vertices'][:]
        observed_slip = h5['vertex_fields/slip'][:,:,:].squeeze()
        h5.close()

        # Sort by y-coordinate.
        reorder = numpy.argsort(observed_coords[:,1])
        self.observed_coords = observed_coords[reorder,:]
        self.observed_slip = observed_slip[reorder,:]

    def _load_raw(self, filename: str):
        data = numpy.loadtxt(filename)
        self.raw_y = data[:,0]
        self.raw_slip = data[:,1]

    def _load_models(self, filename: str):
        model = numpy.loadtxt(filename)
        self.models_coords = model[:,0]
        self.models_slip = model[:, 1:]

        with open(filename, "r") as fin:
            header = fin.readline()
        self.models_penalty = list(map(float, re.findall(r"penalty=(\d+\.\d+)", header)))

    def _plot(self, prescribed_only=False, show_plot=False):
        figure = pyplot.figure(figsize=(6.5, 4.0), layout="tight")
        axes = self._setup_axes(figure)
        self._plot_observed(axes)
        if  not prescribed_only:
            self._plot_models(axes)
            axes.legend(loc="upper right")

        if show_plot:
            pyplot.show()
        filename = "step04-slip.pdf" if prescribed_only else "step06_inversion-results.pdf" 
        figure.savefig(filename)

    def _setup_axes(self, figure):
        axes = figure.add_subplot()
        axes.spines[["top", "right"]].set_visible(False)
        axes.set_xlabel("Distance along strike, km")
        axes.set_ylabel("Left-lateral slip, m")
        axes.set_xlim(self.XLIM)
        return axes

    def _plot_observed(self, axes):
        axes.plot(self.raw_y / km.value, self.raw_slip, linewidth=2, color=self.COLOR_RAW, label="Exact slip")
        axes.plot(self.observed_coords[:,1] / km.value, self.observed_slip[:,1], linewidth=2, color=self.COLOR_OBSERVED, label="Prescribed slip")

    def _plot_models(self, axes):
        nmodels = self.models_slip.shape[1]
        coords = self.models_coords / km.value
        for imodel in range(nmodels):
            axes.plot(coords, self.models_slip[:,imodel], "--", color=self.COLORS_MODEL[imodel], label=f"Model penalty={self.models_penalty[imodel]}")


if __name__ == "__main__":
    cli()


# End of file
