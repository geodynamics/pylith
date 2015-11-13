import h5py
import numpy

def getData(sim):

    filename = "output/%s-domain.h5" % sim
    h5 = h5py.File(filename, "r")
    vertices = h5['geometry/vertices'][:]
    disp = h5['vertex_fields/displacement'][:]
    vel = h5['vertex_fields/velocity'][:]
    timePt = h5['time'][:,0,0]
    h5.close()

    target = (0.0e+3, 1.0e+3)
    dist = (vertices[:,0]-target[0])**2 + (vertices[:,1]-target[1])**2
    ipt = numpy.argmin(dist)
    dispPt = disp[:,ipt,:]
    velPt = vel[:,ipt,:]

    filename = "output/%s-statevars.h5" % sim
    h5 = h5py.File(filename, "r")
    vertices = h5['geometry/vertices'][:]
    cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
    stress = h5['cell_fields/stress'][:]
    stressCauchy = h5['cell_fields/cauchy_stress'][:]
    timeCell = h5['time'][:,0,0]
    h5.close()

    cellCoords = vertices[cells,:]
    cellCenters = numpy.mean(cellCoords, axis=1)
    target = (-500.0, 500.0)
    dist = (cellCenters[:,0]-target[0])**2 + (cellCenters[:,1]-target[1])**2
    icell = numpy.argmin(dist)

    from pyre.units.time import year
    timePt /= year.value
    timeCell /= year.value
    
    return {'time_pt': timePt,
            'disp': dispPt,
            'vel': velPt, 
            'time_cell': timeCell,
            'stress': stress[:,icell,:],
            'cauchy_stress': stressCauchy[:,icell,:]}

# ======================================================================
prob = "extension"

length = 2000.0
width = 0.5*length
p_vs = 3000.0
p_vp = 5196.152422706632
p_density = 2500.0

p_mu = p_density*p_vs**2
p_lambda = p_density*p_vp**2 - 2*p_mu

x0 = 0.0
y0 = -width
x = 0.0
y = 0.0

if prob == "compression":
    vx = -1.0
else:
    vx = +1.0
t = numpy.arange(0.0, 999.01, 1.0)
ex = vx*t / (0.5*length)
ux = ex*(x-x0)
Exx = ex+0.5*ex**2
Eyy = -p_lambda/(p_lambda+2*p_mu)*Exx
ey = -1 + (1 + 2*Eyy)**0.5
uy = ey*(y-y0)
Sxx = p_lambda*(Exx+Eyy) + 2.0*p_mu*Exx
Syy = p_lambda*(Exx+Eyy) + 2.0*p_mu*Eyy

# Elastic
sim = "%s" % prob
dataS = getData(sim)


import pylab
pylab.subplot(2,2,1)
pylab.plot(t, uy, 'k-',
           dataS['time_pt'], dataS['disp'][:,1], 'r--')
pylab.xlabel("Time (year)")
pylab.xlabel("Displacement (m)")
pylab.title("Displacement at top center")

pylab.subplot(2,2,2)
pylab.plot(dataS['time_pt'], dataS['vel'][:,1], 'r--')
pylab.xlabel("Time (year)")
pylab.xlabel("Velocity (m/s)")
pylab.title("Velocity at top center")

pylab.subplot(2,2,3)
pylab.plot(t, Sxx, 'k-',
           dataS['time_cell'], dataS['stress'][:,0], 'r--', 
           dataS['time_cell'], dataS['cauchy_stress'][:,0], 'g:')
pylab.xlabel("Time (year)")
pylab.ylabel("Stress_xx (Pa)")
pylab.title("Stress at mid-depth at left quarter point")

pylab.subplot(2,2,4)
pylab.plot(t, Syy, 'k-',
           dataS['time_cell'], dataS['stress'][:,1], 'r--', 
           dataS['time_cell'], dataS['cauchy_stress'][:,1], 'g:')
pylab.xlabel("Time (year)")
pylab.ylabel("Stress_yy (Pa)")
pylab.title("Stress at mid-depth at left quarter point")

pylab.show()

