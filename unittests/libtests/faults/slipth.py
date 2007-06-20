from math import *

finalSlip = (2.5,  0.0, 0.3)

slipMag = (finalSlip[0]**2+finalSlip[1]**2+finalSlip[2]**2)**0.5
peakRate = 1.6
slipTime = 1.4
t = 2.134

tau = slipMag / (exp(1.0) * peakRate)
t0 = slipTime
slipN = 1.0 - exp(-(t-t0)/tau) * (1.0 + (t-t0)/tau)

print "%13.11f, %13.11f, %13.11f" % \
      (slipN*finalSlip[0], slipN*finalSlip[1], slipN*finalSlip[2])
