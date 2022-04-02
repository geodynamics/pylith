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

p_density = 2500.0
p_vs = 3464.1016
p_vp = 6000.0
p_viscosity = 9.467279917257993e17
p_power_law_exponent = 1.0
p_power_law_reference_strain_rate = 1.0e-6
p_power_law_reference_stress = 2.0*p_viscosity*p_power_law_reference_strain_rate
p_mu = p_density*p_vs*p_vs
p_lambda = p_density*p_vp*p_vp - 2.0*p_mu
p_youngs = p_mu*(3.0*p_lambda + 2.0*p_mu)/(p_lambda + p_mu)
p_poissons = 0.5*p_lambda/(p_lambda + p_mu)

# Uniform stress field (plane strain).
p_t0 = -1.0e+9

params = {
    'p_density': p_density,
    'p_vs': p_vs,
    'p_vp': p_vp,
    'p_viscosity': p_viscosity,
    'p_power_law_exponent': p_power_law_exponent,
    'p_power_law_reference_strain_rate': p_power_law_reference_strain_rate,
    'p_power_law_reference_stress': p_power_law_reference_stress,
    'p_mu': p_mu,
    'p_lambda': p_lambda,
    'p_youngs': p_youngs,
    'p_poissons': p_poissons,
    'T0': p_t0,
    }


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
        gendb.run(params)

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
            if (sim == 'axialtraction_powerlaw_n1'):
                plotter.run(self.axes, isim, params, dt, test_jacobian)
            else:
                plotter.run(self.axes, isim, dt, test_jacobian)

    def _show_plot(self):
        pyplot.show()

    def _parse_command_line(self):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--sim", action="store", dest="sim", default="all", choices=list(self.SIMS.keys()), help="Simulations to run")
        parser.add_argument("--test_jacobian", action="store_true", dest="test_jacobian", default=False, help="Test PyLith Jacobian")

        return parser.parse_args()


# ======================================================================
if __name__ == "__main__":
    RunnerApp().main()


# End of file
