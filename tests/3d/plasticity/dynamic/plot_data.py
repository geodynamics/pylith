#!/usr/bin/env python
import sys

if len(sys.argv) != 3:
    raise ValueError("usage: plot_data.py CELL SIM")

cell = sys.argv[1]
sim = sys.argv[2]

if not sim in ['swave', 'pwave']:
    raise ValueError("Unknown sim '%s'." % sim)


# ======================================================================
import tables
import pylab
import numpy

# ----------------------------------------------------------------------
def swave_soln():

    dt = 5.0e-6
    tend = 5.0e-4
    t = numpy.arange(0.0, tend+0.5*dt, dt)

    mask = t < 1.2956829537e-4
    sxx = 0*t
    syy = 0*t
    szz = 0*t
    sxy = mask*2.9403000000e+10*t + ~mask*3.8096965888e+6
    syz = 0*t
    sxz = 0*t

    return {'t': t,
            'sxx': sxx,
            'syy': syy,
            'szz': szz,
            'sxy': sxy,
            'syz': syz,
            'sxz': sxz}


# ----------------------------------------------------------------------
def pwave_soln():
    dt = 5.0e-5
    tend = 5.0e-3
    t = numpy.arange(0.0, tend+0.5*dt, dt)

    mask = t < 1.7246213245e-3
    sxx = mask*-8.8216171200e+10*t + ~mask*(-8.5665432299e+10*t-4.3990587021e+6)
    syy = mask*-2.9410171200e+10*t + ~mask*(-3.0685540651e+10*t+2.1995293511e+6)
    szz = 1.0*syy
    sxy = 0*t
    syz = 0*t
    sxz = 0*t

    return {'t': t,
            'sxx': sxx,
            'syy': syy,
            'szz': szz,
            'sxy': sxy,
            'syz': syz,
            'sxz': sxz}


# ----------------------------------------------------------------------
# Analytical soln
if sim == "swave":
    analytic = swave_soln()
elif sim == "pwave":
    analytic = pwave_soln()
else:
    raise ValueError("Unknown sim.")


# Synthetic
h5 = tables.openFile("%s_%s-statevars.h5" % (cell, sim), "r")
time = h5.root.time[:].ravel()
stress = h5.root.cell_fields.stress[:]
h5.close()

# Redimensionalize
timescale = 1.0e-6
analytic['t'] /= timescale
time /= timescale

stressscale = 1.0e+6
analytic['sxx'] /= stressscale
analytic['syy'] /= stressscale
analytic['szz'] /= stressscale
analytic['sxy'] /= stressscale
analytic['syz'] /= stressscale
analytic['sxz'] /= stressscale
stress /= stressscale

# Plot
fig = pylab.Figure()

if sim == "swave":
    xlim = (0, 500)
    ylim = (-1, 10)
else:
    xlim = (0, 5000)
    ylim = (-450, 10)


# Axial components
ax = pylab.subplot(1, 2, 1)
ax.plot(analytic['t'], analytic['sxx'], 'k-',
        time, stress[:,0,0:48:6], 'r--',
        analytic['t'], analytic['syy'], 'k-',
        time, stress[:,0,1:48:6], 'r--',
        analytic['t'], analytic['szz'], 'k-',
        time, stress[:,0,2:48:6], 'r--')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xlabel("Time (usec)")
ax.set_ylabel("Stress (MPa)")

# Shear components
ax = pylab.subplot(1, 2, 2)
ax.plot(analytic['t'], analytic['sxy'], 'k-',
        time, stress[:,0,3:48:6], 'r--',
        analytic['t'], analytic['syz'], 'k-',
        time, stress[:,0,4:48:6], 'r--',
        analytic['t'], analytic['sxz'], 'k-',
        time, stress[:,0,5:48:6], 'r--')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xlabel("Time (usec)")
ax.set_ylabel("Stress (MPa)")

pylab.show()
