#!/usr/bin/env python

sim = "ratestate_stable"

# ======================================================================
import tables
import pylab
import numpy
from math import exp

# ----------------------------------------------------------------------
dt = 0.01
t = numpy.arange(0.0, 12.001, dt)
mu0 = 0.6
a = 0.016
b = 0.012
L = 2.0e-6
V0 = 1.0e-6

def integrateStateVar(theta, v, t0):
    pts = numpy.where(t > t0)[0]
    for i in pts:
        theta[i] = theta[i-1]*exp(-v/L*dt) + \
            dt - 0.5*(v/L)*dt**2 + 1.0/6.0*(v/L)*dt**3;
    return


theta = L/V0*numpy.ones(t.shape)

mask1 = t < 2.0
V1 = 1.0e-6
integrateStateVar(theta, V1, 0.0)

mask2 = numpy.bitwise_and(t >= 2.0, t < 4.0)
V2 = 2.0e-5
integrateStateVar(theta, V2, 2.0)

mask3 = numpy.bitwise_and(t >= 4.0, t < 8.0)
V3 = 5.0e-6
integrateStateVar(theta, V3, 4.0)

mask4 = numpy.bitwise_and(t >= 8.0, t < 12.0)
V4 = 2.0e-7
integrateStateVar(theta, V4, 8.0)

slipRateE = mask1*V1 + mask2*V2 + mask3*V3 + mask4*V4
stateVarE = theta

muE = mu0 + a*numpy.log(slipRateE/V0) + b*numpy.log(V0*stateVarE/L)

# ----------------------------------------------------------------------

h5 = tables.openFile("output/%s-fault.h5" % sim, "r")
time = h5.root.time[:].ravel()
slip = h5.root.vertex_fields.slip[:]
slipRate = h5.root.vertex_fields.slip_rate[:]
stateVar = h5.root.vertex_fields.state_variable[:]
traction = h5.root.vertex_fields.traction[:]
h5.close()

fig = pylab.Figure()

p = 2

ax = pylab.subplot(1, 4, 1)
ax.plot(time, slip[:,p,0])

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
