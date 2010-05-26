from math import *

finalSlip = (2.6, -1.0, 0.4)
riseTime = 1.8

slipMag = (finalSlip[0]**2+finalSlip[1]**2+finalSlip[2]**2)**0.5

slipTime1 = 2.0
slipTime2 = 0.7
t = 2.134
dt = 0.01

peakRate = slipMag / riseTime * 1.745;
tau = slipMag / (exp(1.0) * peakRate)

t0 = slipTime1
slipN1 = 1.0 - exp(-(t-t0)/tau) * (1.0 + (t-t0)/tau)
slipNp = 1.0 - exp(-(t-dt-t0)/tau) * (1.0 + (t-dt-t0)/tau)
slipNincr1 = slipN1 - slipNp

t0 = slipTime2
slipN2 = 1.0 - exp(-(t-t0)/tau) * (1.0 + (t-t0)/tau)
slipNp = 1.0 - exp(-(t-dt-t0)/tau) * (1.0 + (t-dt-t0)/tau)
slipNincr2 = slipN2 - slipNp

slipN = slipN1 + slipN2
slipNincr = slipNincr1 + slipNincr2

print "%13.11f, %13.11f, %13.11f" % \
      (slipN1*finalSlip[0], slipN1*finalSlip[1], slipN1*finalSlip[2])

print "%13.11f, %13.11f, %13.11f" % \
      (slipNincr*finalSlip[0], slipNincr*finalSlip[1], slipNincr*finalSlip[2])
