import h5py
import numpy

def getData(sim):

    filename = "%s-domain.h5" % sim
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

    filename = "%s-statevars.h5" % sim
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
    print icell

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
prob = "compression"

# Elastic
sim = "%s_dx250" % prob
dataF = getData(sim)

# Viscoelastic
sim = "%s_dx200" % prob
dataI = getData(sim)

import pylab
pylab.subplot(2,2,1)
pylab.plot(dataF['time_pt'], dataF['disp'][:,1], 'r-', 
           dataI['time_pt'], dataI['disp'][:,1], 'b--')
pylab.xlabel("Time (year)")
pylab.xlabel("Displacement (m)")
pylab.title("Displacement at top center")

pylab.subplot(2,2,2)
pylab.plot(dataF['time_pt'], dataF['vel'][:,1], 'r-', 
           dataI['time_pt'], dataI['vel'][:,1], 'b--')
pylab.xlabel("Time (year)")
pylab.xlabel("Velocity (m/s)")
pylab.title("Velocity at top center")

pylab.subplot(2,2,3)
pylab.plot(dataF['time_cell'], dataF['stress'][:,0], 'r:', 
           dataF['time_cell'], dataF['cauchy_stress'][:,0], 'g-', 
           dataI['time_cell'], dataI['stress'][:,0], 'b:',
           dataI['time_cell'], dataI['cauchy_stress'][:,0], 'm--')
pylab.xlabel("Time (year)")
pylab.ylabel("Stress_xx (Pa)")
pylab.title("Stress at mid-depth at left quarter point")

pylab.subplot(2,2,4)
pylab.plot(dataF['time_cell'], dataF['stress'][:,1], 'r:', 
           dataF['time_cell'], dataF['cauchy_stress'][:,1], 'g-', 
           dataI['time_cell'], dataI['stress'][:,1], 'b:',
           dataI['time_cell'], dataI['cauchy_stress'][:,1], 'm--')
pylab.xlabel("Time (year)")
pylab.ylabel("Stress_yy (Pa)")
pylab.title("Stress at mid-depth at left quarter point")

pylab.show()

