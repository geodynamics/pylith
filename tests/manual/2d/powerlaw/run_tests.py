#!/usr/bin/env python

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
#
# @file tests/manual/2d/powerlaw/run_tests.py
#
# @brief Run power-law manual tests.
#

import numpy
import subprocess
import importlib
import matplotlib.pyplot as plt

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
T0 = -1.0e9

params = {'p_density': p_density,
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
          'T0': T0}

allSims = ['axialtraction_powerlaw', 'axialtraction_powerlaw_n1', 'sheartraction_powerlaw']
testJacobianArgs = ['--petsc.snes_test_jacobian', '--petsc.snes_test_jacobian_view']
dt = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]
numSteps = len(dt)
simArgs = {'axialtraction_powerlaw': ['pylith', 'axialtraction_powerlaw.cfg'],
           'axialtraction_powerlaw_n1': ['pylith', 'axialtraction_powerlaw_n1.cfg'],
           'sheartraction_powerlaw': ['pylith', 'sheartraction_powerlaw.cfg']}

def run_sims(sims, jacobianArgs):
    """
    Run simulations and plot results for power-law.
    """

    # Generate spatialdb for axial traction (n=1).
    gendb.run(params)

    # Create plot.
    numSims = len(sims)
    numRows = 4
    useJacobian = True
    if len(jacobianArgs) == 0:
        numRows = 3
        useJacobian = False
    fig, axs = plt.subplots(numRows, numSims)
    axs = axs.reshape(numRows, numSims)
    simNum = 0

    for sim in sims:
        for stepNum in range(numSteps):
            dtStr = repr(dt[stepNum])
            baseName = sim + '_dt' + dtStr
            logFile = baseName + '.log'
            nameArg = '--problem.defaults.name=' + baseName
            dumpArg = '--dump_parameters.filename=output/' + baseName + '-parameters.json'
            progArg = '--problem.progress_monitor.filename=output/' + baseName + '-progress.txt'
            stepArg = '--problem.initial_dt=' + dtStr + '*year'
            metaArg = '--metadata.arguments=[' + simArgs[sim][-1] + ']'
            runargs = simArgs[sim] + [nameArg, dumpArg, progArg, stepArg, metaArg] + jacobianArgs
            f = open(logFile, 'w')
            subprocess.run(runargs, stdout=f, stderr=f, check=True)
            f.close()
        plotModule = 'plot_' + sim
        plot = importlib.import_module(plotModule)
        if (sim == 'axialtraction_powerlaw_n1'):
            plot.run(axs, simNum, params, dt, useJacobian)
        else:
            plot.run(axs, simNum, dt, useJacobian)

        simNum += 1

    plt.show()

# ======================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sims", action="store", dest="sims", default="all", choices=allSims, help="simulations to run")
    parser.add_argument("-j", "--test_jacobians", action="store_true", dest="test_jacobians", default=False, help="test PyLith Jacobians")

    args = parser.parse_args()
    if args.sims == 'all':
        sims = allSims
    else:
        sims = [args.sims]

    jacobianArgs = []
    if args.test_jacobians:
        jacobianArgs = testJacobianArgs

    run_sims(sims, jacobianArgs)
        
# End of file
