# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""Application for visualizing PyLith output."""

import argparse
import pathlib
import h5py

import numpy
import pyvista
from pylith.viz import (PlotMesh, PlotField, WarpGrid, io)


class App:

    def run(self, args):
        pyvista.set_plot_theme("dark")
        filenames = args.filenames.split(",")
        self._verify_filenames(filenames)
        self.data = [io.PyLithReader.read(filename) for filename in filenames]
        if not self.data:
            return

        self.plotter = pyvista.Plotter()
        self.figure = args.func(plotter=self.plotter, data=self.data, args=args)
        
        times = self.data[0].time
        if len(times) > 1 and self.figure.TIME_DEPENDENT:
            self.t_index = 0
            self._add_timestamp_widget(timestamp=times[self.t_index])
            self._add_time_slider()
        self._add_keybindings()

        self.plotter.show_axes()
        self.plotter.view_xy()
        self.plotter.show()

    def _verify_filenames(self, filenames: list):
        missing_files = [filename for filename in filenames if not pathlib.Path(filename).is_file()]
        if len(missing_files):
            raise IOError(f"Could not find PyLith HDF5 output files {missing_files}.")
        invalid_files = []
        for filename in filenames:
            try:
                h5py.File(filename, "r")
            except IOError:
                invalid_files.append(filename)
        if len(invalid_files):
            raise IOError(f"Could not open PyLith HDF5 output files {invalid_files}.")

    def _get_time_label(self, timestamp):
        YEAR = 365.25*24*3600
        return f"Time={timestamp:7.2f} s" if timestamp < 1.0e+5 else f"Time={timestamp / YEAR:7.2f} yr"

    def _get_timestamp(self, t):
        times = self.data[0].time
        return numpy.asarray(times >= t).nonzero()[0][0] if t < times.max() else len(times)-1

    def _add_timestamp_widget(self, timestamp: float):
        t_label = self._get_time_label(timestamp)
        self.timestamp_widget = self.plotter.add_text(t_label, position="upper_left", font_size=12, name="timestamp")

    def _add_time_slider(self):
        times = self.data[0].time
        t_range = (times.min(), times.max())
        self.plotter.add_slider_widget(lambda t: self._update_time(self._get_timestamp(t)), rng=t_range, pointa=(0.35, 0.92), pointb=(0.65, 0.92), value=t_range[0], title="Time")

    def _update_time(self, t_index):
        times = self.data[0].time
        self.t_index = max(0, min(t_index, len(times)-1))
        # Update time label _before_ updating figure to ensure updated label is rendered.
        t_label = self._get_time_label(times[self.t_index])
        self.timestamp_widget.set_text("upper_left", t_label)
        self.figure.update_time(t_index=self.t_index)

    def _add_keybindings(self):
        def increment_time():
            self._update_time(self.t_index+1)
        def decrement_time():
            self._update_time(self.t_index-1)

        lines = []
        times = self.data[0].time
        if len(times) > 1 and self.figure.TIME_DEPENDENT:
            self.plotter.add_key_event("p", decrement_time)
            self.plotter.add_key_event("n", increment_time)
            lines += ["    'p': Decrement time"]
            lines += ["    'n': Increment time"]
        if len(lines):
            print("Key bindings:")
            print("\n".join(lines))


def cli():
    """Parse command line arguments and run application.
    """

    parser = argparse.ArgumentParser(
        description="Application for visualizing PyLith output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--filenames",
        action="store",
        required=True,
        help="Comma separated list of PyLith HDF output files."
    )
    subparsers = parser.add_subparsers(title="subcommands", required=True, dest="subcommand")

    PlotMesh.add_args(subparsers.add_parser("plot_mesh"))
    WarpGrid.add_args(subparsers.add_parser("warp_grid"))
    PlotField.add_args(subparsers.add_parser("plot_field"))

    args = parser.parse_args()
    App().run(args)


# End of file
