#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    raise ValueError("usage: plot_friction.py [weak | stable]")

sim = sys.argv[1]

if not sim in ['weak', 'stable']:
    raise ValueError("Unknown sim '%s'." % sim)

# ======================================================================
import h5py
import pylab
import numpy
from math import exp

# ----------------------------------------------------------------------
dt = 0.01
t = numpy.arange(0.0, 14.001, dt)
mu0 = 0.6
if sim == "stable":
    a = 0.016
    b = 0.012
elif sim == "weak":
    a = 0.008
    b = 0.012
else:
    raise ValueError("Unknown sim '%s'." % sim)
L = 2.0e-6
V0 = 1.0e-6

def integrateStateVar(theta, v, t0):
    pts = numpy.where(t > t0)[0]
    for i in pts:
        theta[i] = theta[i-1]*exp(-v/L*dt) + \
            dt - 0.5*(v/L)*dt**2 + 1.0/6.0*(v/L)*dt**3;
    return



mask1 = t < 2.0
V1 = 1.0e-6
theta = L/V1*numpy.ones(t.shape)
integrateStateVar(theta, V1, 0.0)

mask2 = numpy.bitwise_and(t >= 2.0, t < 4.0)
V2 = 1.0e-5
integrateStateVar(theta, V2, 2.0)

mask3 = numpy.bitwise_and(t >= 4.0, t < 6.0)
V3 = 4.0e-6
integrateStateVar(theta, V3, 4.0)

mask4 = numpy.bitwise_and(t >= 6.0, t < 8.0)
V4 = 2.0e-5
integrateStateVar(theta, V4, 6.0)

mask5 = numpy.bitwise_and(t >= 8.0, t < 12.0)
V5 = 5.0e-6
integrateStateVar(theta, V5, 8.0)

mask6 = t >= 12.0
V6 = 1.0e-6
integrateStateVar(theta, V6, 12.0)

slipRateE = mask1*V1 + mask2*V2 + mask3*V3 + mask4*V4 + mask5*V5 + mask6*V6
stateVarE = theta

muE = mu0 + a*numpy.log(slipRateE/V0) + b*numpy.log(V0*stateVarE/L)

# ----------------------------------------------------------------------

h5 = h5py.File("output/ratestate_%s-fault.h5" % sim, "r")
time = h5['time'][:].ravel()
slip = h5['vertex_fields/slip'][:]
slipRate = h5['vertex_fields/slip_rate'][:]
stateVar = h5['vertex_fields/state_variable'][:]
traction = h5['vertex_fields/traction'][:]
h5.close()

fig = pylab.Figure()

p = 2

ax = pylab.subplot(1, 4, 1)
ax.plot(time, slip[:,p,0], 'r--')

ax = pylab.subplot(1, 4, 2)
ax.plot(t, numpy.log10(numpy.abs(slipRateE)), 'b-',
        time, numpy.log10(numpy.abs(slipRate[:,p,0])), 'r--')
ax.set_ylim(-12, 0.0)

ax = pylab.subplot(1, 4, 3)
ax.plot(t, numpy.log10(stateVarE), 'b-',
        time, numpy.log10(stateVar[:,p,0]), 'r--')
ax.set_ylim(-2.0,6.0)

ax = pylab.subplot(1, 4, 4)
ax.plot(t, muE, 'b-',
        time, numpy.fabs(traction[:,p,0]/traction[:,p,1]), 'r--')


pylab.show()
