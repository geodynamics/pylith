#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

sim = 'genmaxwell_QpQs'
problem = 'bulk'
dt = 1.0

# ======================================================================
import tables
import numpy
import pylab
from mypylab.Figure import Figure


vp = 5291.502622129181
vs = 3000.0
density = 2500.0
mu0 = density*vs**2
k0 = density*vp**2 - 4*mu0/3.0
strain0 = numpy.array([0.0, 0.0, +1.0/10.0e+3, 0.0, 0.0, 0.0],
                      dtype=numpy.float64)

if problem == 'shear':
    tm_shear = [25, 50, 100]
    tm_bulk = [25, 50, 100]
    ratio_shear = [0.5, 0.2, 0.0]
    ratio_shear0 = 1.0 - sum(ratio_shear)
    ratio_bulk = [0.0, 0.0, 0.0]
    ratio_bulk0 = 1.0 - sum(ratio_bulk)
elif problem == 'bulk':
    tm_shear = [25, 50, 100]
    tm_bulk = [20, 40, 80]
    ratio_shear = [0.4, 0.1, 0.2]
    ratio_shear0 = 1.0 - sum(ratio_shear)
    ratio_bulk = [0.1, 0.2, 0.3]
    ratio_bulk0 = 1.0 - sum(ratio_bulk)
    
# ----------------------------------------------------------------------
filename = "output/%s_%s-statevars.h5" % (problem, sim)

h5 = tables.openFile(filename, 'r')
stress = h5.root.cell_fields.stress[:] / 1.0e+6
(ntimesteps, npts, tensorSize) = stress.shape
h5.close()    

#time =  h5.root.vertex_fields.time (not yet available)
t = numpy.arange(0, dt*ntimesteps, dt)
# END TEMPORARY

mu_t = mu0 * (ratio_shear0 + \
                  ratio_shear[0]*numpy.exp(-t/tm_shear[0]) + \
                  ratio_shear[1]*numpy.exp(-t/tm_shear[1]) + \
                  ratio_shear[2]*numpy.exp(-t/tm_shear[2]))
k_t = k0 * (ratio_bulk0 + \
                ratio_bulk[0]*numpy.exp(-t/tm_bulk[0]) + \
                ratio_bulk[1]*numpy.exp(-t/tm_bulk[1]) + \
                ratio_bulk[2]*numpy.exp(-t/tm_bulk[2]))
lambda_t = k_t - 2.0/3.0*mu_t
stressE = numpy.zeros( (ntimesteps,6), dtype=numpy.float64)
stressE[:,0] = lambda_t*numpy.sum(strain0[0:3]) + 2.0*mu_t*strain0[0]
stressE[:,1] = lambda_t*numpy.sum(strain0[0:3]) + 2.0*mu_t*strain0[1]
stressE[:,2] = lambda_t*numpy.sum(strain0[0:3]) + 2.0*mu_t*strain0[2]
stressE[:,3] = 2.0*mu_t*strain0[3]
stressE[:,4] = 2.0*mu_t*strain0[4]
stressE[:,5] = 2.0*mu_t*strain0[5]
stressE /= 1.0e+6

print stress[0:2,0,:]
print stress[ntimesteps-1,0,:]
print stressE[0,:]
print stressE[1,:]
print stressE[ntimesteps-1,:]

# ----------------------------------------------------------------------
nrows = 2
ncols = 3
irow = 1
icol = 1

fig = Figure(fontsize=8, color="lightbg")
fig.open(9.0, 7.0, margins=[[0.45, 0.25, 0.1],
                            [0.45, 0.5, 0.2]])


icomp = 0
for irow in xrange(1, 3):
  for icol in xrange(1, 4):
      ax = fig.axes(nrows, ncols, irow, icol)
      
      ax.plot(t, stressE[:,icomp], 'r-',
              t, stress[:,0,icomp], 'b--')
      #ax.plot(t, devStress[:,0,icomp])
      ax.set_ylim( (-2, 8) )

      icomp += 1

pylab.show()


# End of file
