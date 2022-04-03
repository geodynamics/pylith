#!/usr/bin/env nemesis

# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2018 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

import contextlib
import io
import importlib
import matplotlib.pyplot as pyplot

from pylith.apps.PyLithApp import PyLithApp

import axialtraction_powerlaw_n1_gendb as gendb


class RunnerApp():

    SIMS = {
        'axialtraction_powerlaw': ['axialtraction_powerlaw.cfg'],
        'axialtraction_powerlaw_n1': ['axialtraction_powerlaw_n1.cfg'],
        'sheartraction_powerlaw': ['sheartraction_powerlaw.cfg'],
    }
    TIMESTEPS = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]

    def main(self):
        """Application entry point.
        """
        args = self._parse_command_line()
        if args.sim == 'all':
            sims = self.SIMS
        else:
            sims = {args.sim: self.SIMS[args.sim]}

        self._generate_db()
        self._setup_plot(len(sims), args.test_jacobian)
        self._run_sims(sims, args.test_jacobian)
        self._show_plot()

    def _generate_db(self):
        gendb.run()

    def _setup_plot(self, nsims, use_jacobian):
        nrows = 4 if use_jacobian else 3
        fig, axes = pyplot.subplots(nrows, nsims)
        self.axes = axes.reshape(nrows, nsims)

    def _run_sims(self, sims, test_jacobian):
        ARGS_JACOBIAN = ['--petsc.snes_test_jacobian', '--petsc.snes_test_jacobian_view']
        args_jacobian = ARGS_JACOBIAN if test_jacobian else []

        for isim, sim in enumerate(sims.keys()):
            for dt in self.TIMESTEPS:
                sim_fullname = f"{sim}_dt{dt:04.2f}"
                log_filename = f"{sim_fullname}.log"

                args_sim = [
                    f"--problem.defaults.name={sim_fullname}",
                    f"--dump_parameters.filename=output/{sim_fullname}-parameters.json",
                    f"--problem.progress_monitor.filename=output/{sim_fullname}-progress.txt",
                    f"--problem.initial_dt={dt:04.2f}*year",
                ]
                pylith_args = sims[sim] + args_sim + args_jacobian
                with contextlib.redirect_stdout(io.StringIO()) as sout:
                    app = PyLithApp()
                    app.run(argv=["pylith"] + pylith_args)
                with open(log_filename, "w") as log:
                    log.write(sout.getvalue())

            plot_module = f"plot_{sim}"
            plotter = importlib.import_module(plot_module)
            plotter.run(self.axes, isim, dt, test_jacobian)

    def _show_plot(self):
        pyplot.show()

    def _parse_command_line(self):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--sim", action="store", dest="sim", default="all", choices=list(self.SIMS.keys()), help="Simulations to run")
        parser.add_argument("--test_jacobian", action="store_true", dest="test_jacobian", default=False, help="Test PyLith Jacobian")

        return parser.parse_args()


if __name__ == "__main__":
    RunnerApp().main()


# End of file
