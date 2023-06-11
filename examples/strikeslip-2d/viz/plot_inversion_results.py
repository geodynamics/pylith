#!/usr/bin/env nemesis
"""
This script generates a plot comparing the predicted solutions to
the true solution for Step 6.
"""

# Import argparse and re (regular expressions) Python modules(standard Python)
import argparse
import re

# Import numpy, h5py, and matplotlib modules (included in PyLith binary installation)
import numpy
import h5py
import matplotlib.pyplot as pyplot

class PlotApp:

    FILENAME_MODELS = "output/step06_inversion-results.txt"
    FILENAME_OBSERVED = "output/step04_varslip-fault.h5"

    
    COLOR_OBSERVED = "black"
    COLORS_MODEL = ["blue", "orange", "purple", "red", "cyan"]
    XLIM = (-25.0, 25.0)

    def main(self, filename_models: str=None, filename_observed: str=None, observed_only: bool=False, show_plot: bool=False):
        filename_models = filename_models or self.FILENAME_MODELS
        filename_observed = filename_observed or self.FILENAME_OBSERVED

        self.load_observed(filename_observed)
        self.load_models(filename_models)
        self.plot(observed_only, show_plot)

    def load_observed(self, filename):
        h5 = h5py.File(filename, "r")
        observed_coords = h5['geometry/vertices'][:]
        observed_slip = h5['vertex_fields/slip'][:,:,:].squeeze()
        h5.close()

        # Sort by y-coordinate.
        reorder = numpy.argsort(observed_coords[:,1])
        self.observed_coords = observed_coords[reorder,:]
        self.observed_slip = observed_slip[reorder,:]

    def load_models(self, filename):
        model = numpy.loadtxt(filename)
        self.models_coords = model[:,0]
        self.models_slip = model[:, 1:]

        with open(filename, "r") as fin:
            header = fin.readline()
        self.models_penalty = list(map(float, re.findall(r"penalty=(\d+\.\d+)", header)))

    def plot(self, observed_only=False, show_plot=False):
        figure = pyplot.figure(figsize=(7.0, 3.5), dpi=150, layout="tight")
        axes = self._setup_axes(figure)
        self._plot_observed(axes)
        if  not observed_only:
            self._plot_models(axes)
            axes.legend(loc="upper right")

        if show_plot:
            pyplot.show()
        filename = "step04-slip.pdf" if observed_only else "step06_inversion-results.pdf" 
        figure.savefig(filename)

    def _setup_axes(self, figure):
        axes = figure.add_subplot()
        axes.set_xlabel("Distance along Strike (km)")
        axes.set_ylabel("Left-lateral Slip (m)")
        axes.set_xlim(self.XLIM)
        return axes

    def _plot_observed(self, axes):
        axes.plot(self.observed_coords[:,1] / 1.0e+3, self.observed_slip[:,1], linewidth=2, color=self.COLOR_OBSERVED, label="Observed")

    def _plot_models(self, axes):
        nmodels = self.models_slip.shape[1]
        coords = self.models_coords / 1.0e+3
        for imodel in range(nmodels):
            axes.plot(coords, self.models_slip[:,imodel], color=self.COLORS_MODEL[imodel], label=f"Model penalty={self.models_penalty[imodel]}")


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", action="store", dest="filename_models", type=str, default=PlotApp.FILENAME_MODELS, help="Filename of output from inversion.")
    parser.add_argument("--observed", action="store", dest="filename_observed", type=str, default=PlotApp.FILENAME_OBSERVED, help="Name of HDF5 file with station output with fake observations.")
    parser.add_argument("--observed-only", action="store_true", dest="observed_only", default=False, help="Show only observed slip.")
    parser.add_argument("--no-gui", action="store_false", dest="show_plot", default=True, help="Do not display plot.")

    args = parser.parse_args()
    kwargs = {
        "filename_models": args.filename_models,
        "filename_observed": args.filename_observed,
        "observed_only": args.observed_only,
        "show_plot": args.show_plot,
        }
    PlotApp().main(**kwargs)
    
    
if __name__ == "__main__":
    cli()


# End of file
