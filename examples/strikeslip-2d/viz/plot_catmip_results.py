#!/usr/bin/env nemesis
"""
This script generates a plot comparing the predicted solutions to
the true solution for Step 7 (CATMIP inversions).
"""

# Standard Python modules
import pathlib

# Python modules included in PyLith binary installation
import numpy
import h5py
import matplotlib.pyplot as pyplot
from pythia.pyre.units.length import km

OUTPUT_DIR = pathlib.Path("output")

FILENAME_IMPULSES = OUTPUT_DIR / "step05_greensfns-fault.h5"
FILENAME_PRESCRIBED = OUTPUT_DIR / "step04_varslip-fault.h5"
FILENAME_RAW = OUTPUT_DIR / "slip_variable.txt"


def cli():
    """Command line interface.
    """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--catmip-theta", action="store", dest="filename_theta", type=str, required=True, help="Filename of output from CATMIP inversion.")
    parser.add_argument("--pylith-impulses", action="store", dest="filename_impulses", type=str, default=FILENAME_IMPULSES, help="Name of HDF5 file with fault slip inpulses from Green's functions simulation.")
    parser.add_argument("--prescribed-slip", action="store", dest="filename_prescribed", type=str, default=FILENAME_PRESCRIBED, help="Name of HDF5 file for fault from prescribed slip simulation.")
    parser.add_argument("--raw", action="store", dest="filename_raw", type=str, default=FILENAME_RAW, help="Name of ASCII file with raw slip profile.")
    parser.add_argument("--no-gui", action="store_false", dest="show_plot", default=True, help="Do not display plot.")

    args = parser.parse_args()
    PlotApp().run(
        filename_prescribed=args.filename_prescribed,
        filename_impulses=args.filename_impulses,
        filename_theta=args.filename_theta,
        filename_raw=args.filename_raw,
        show_plot=args.show_plot,
        )

            
class PlotApp:
    """Application to plot CATMIP inversion results.
    """
    COLOR_RAW = "gray"
    COLOR_OBSERVED = "black"
    COLOR_MODEL = "orange"
    XLIM = (-25.0, 35.0)

    def run(self, filename_theta: str, filename_prescribed: str, filename_impulses: str, filename_raw: str, show_plot: bool):
        """Run plotting application.

        Args:
            filename_theta: Filename of output from inversion.
            filename_prescribed: Name of HDF5 file with station output with fake observations.
            filename_raw: Name of ASCII file with raw slip values.
            prescribed_only: Show only prescribed slip.
            show_plot: Show plot window.
        """
        self._load_observed(filename_prescribed)
        self._load_raw(filename_raw)
        self._load_impulses(filename_impulses)
        self._load_inversion_results(filename_theta)

        sim_stem = pathlib.Path(filename_theta).name.split("theta")[0]
        filename_plot = sim_stem + "results.pdf"
        
        slip_median = self._calc_median()
        slip_stddev = self._calc_stddev()
        slip_minmax = self._calc_minmax()
        self._plot(slip_median, slip_stddev, slip_minmax, filename_plot, show_plot)

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

    def _load_impulses(self, filename: str):
        h5 = h5py.File(filename, "r")
        impulse_coords = h5['geometry/vertices'][:]
        impulse_slip = h5['vertex_fields/slip'][:,:,:].squeeze()
        h5.close()

        # Sort by y-coordinate.
        reorder = numpy.argsort(impulse_coords[:,1])
        self.impulse_coords = impulse_coords[reorder,:]
        self.impulse_slip = impulse_slip[:,reorder,:]

    def _load_inversion_results(self, filename: str):
        with open(filename, "rb") as fin:
            coefs = numpy.frombuffer(fin.read(), dtype=numpy.float64)
        nimpulses = self.impulse_slip.shape[0]
        self.inversion_coefs = coefs.reshape((nimpulses,-1))

    def _calc_median(self) -> numpy.ndarray:
        median_coefs = numpy.median(self.inversion_coefs, axis=1)
        return self._calc_slip(median_coefs)
        
    def _calc_stddev(self) -> numpy.ndarray:
        median_coefs = numpy.median(self.inversion_coefs, axis=1)
        stddev_coefs = numpy.std(self.inversion_coefs, axis=1)
        low_slip = self._calc_slip(median_coefs - stddev_coefs)
        high_slip = self._calc_slip(median_coefs + stddev_coefs)
        return (low_slip, high_slip)
        
    def _calc_minmax(self) -> tuple:
        min_coefs = numpy.min(self.inversion_coefs, axis=1) 
        max_coefs = numpy.max(self.inversion_coefs, axis=1) 
        min_slip = self._calc_slip(min_coefs)
        max_slip = self._calc_slip(max_coefs)
        return (min_slip, max_slip)
        
    def _calc_slip(self, slip_coefs) -> numpy.ndarray:
        nimpulses, nlocs = self.impulse_slip.shape[0:2]
        slip = numpy.dot(slip_coefs.reshape((1,-1)), self.impulse_slip.reshape((nimpulses,-1)))
        return slip.reshape((nlocs,-1))

    def _plot(self, median: numpy.ndarray, stddev: numpy.ndarray, minmax: tuple, filename: str, show_plot: bool):
        figure = pyplot.figure(figsize=(6.5, 4.0), layout="tight")
        axes = self._setup_axes(figure)
        self._plot_observed(axes)
        self._plot_median(median, axes)
        self._plot_stddev(stddev, axes)
        self._plot_minmax(minmax, axes)

        axes.legend(loc="upper right")

        if show_plot:
            pyplot.show()
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
        axes.plot(self.observed_coords[:,1] / 1.0e+3, self.observed_slip[:,1], linewidth=2, color=self.COLOR_OBSERVED, label="Prescribed slip")

    def _plot_median(self, median, axes):
        axes.plot(self.impulse_coords[:,1] / 1.0e+3, median[:,1], linewidth=2, color=self.COLOR_MODEL, label="Model median")

    def _plot_stddev(self, stddev, axes):
        low_slip, high_slip = stddev
        axes.fill_between(self.impulse_coords[:,1] / 1.0e+3, low_slip[:,1], high_slip[:,1], color=self.COLOR_MODEL, lw=0, alpha=0.5, label=r"median$\pm$stddev")

    def _plot_minmax(self, minmax, axes):
        min_slip, max_slip = minmax
        axes.plot(self.impulse_coords[:,1] / 1.0e+3, min_slip[:,1], color=self.COLOR_MODEL, lw=1, ls="--", label="min-max")
        axes.plot(self.impulse_coords[:,1] / 1.0e+3, max_slip[:,1], color=self.COLOR_MODEL, lw=1, ls="--", label=None)


if __name__ == "__main__":
    cli()


# End of file
