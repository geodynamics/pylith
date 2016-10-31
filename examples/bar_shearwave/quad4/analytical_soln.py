#!/usr/bin/env python

# Compute the analytical steady-state solution.
#
# The shear traction on the fault is balanced by the shear traction
# from friction and the traction due to deformation of the
# bar. The absorbing dampers apply a traction on the end of the bar
# equal to the shear wave speed times the density.

# sim = {'static', 'slipweak', 'ratestate'}
sim = "ratestate"

if sim == "static":
    mu_f = 0.6
elif sim == "slipweak":
    mu_f = 0.5
elif sim == "ratestate":
    L = 0.02 # Dc
    mu_f = 0.5313
else:
    raise ValueError("Unknown sim '%s'." % sim)

# Traction applied to the fault
Ts = 75.0e+6 + 25.0e+6

# Shear traction due to friction
Tf = 120.0e+6 * mu_f

# Shear wave speed and density from matprops.spatialdb
Vs = 1000.0
density = 2500.0

# Compute the slip rate on the fault
sliprate = 2.0*(Ts-Tf) / (Vs*density)

# Compute the shear strain in the bar
shearmodulus = density*Vs**2
shearstrain = (Ts-Tf) / (2.0*shearmodulus)

print "Fault traction (friction): ", Tf/1.0e+6, "MPa"
print "Slip rate: ", sliprate, "m/s"
print "Shear strain: ", shearstrain
if sim == "ratestate":
    print "Rate-state friction state var: ", L/sliprate, "s"
