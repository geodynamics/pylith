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

    target = (0.0e+3, 2.0e+3)
    dist = (vertices[:,0]-target[0])**2 + (vertices[:,1]-target[1])**2
    ipt = numpy.argmin(dist)
    dispPt = disp[:,ipt,:]
    velPt = vel[:,ipt,:]

    filename = "%s-elastic.h5" % sim
    h5 = h5py.File(filename, "r")
    vertices = h5['geometry/vertices'][:]
    cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
    stress = h5['cell_fields/stress'][:]
    stressCauchy = h5['cell_fields/cauchy_stress'][:]
    timeCell = h5['time'][:,0,0]
    h5.close()

    cellCoords = vertices[cells,:]
    cellCenters = numpy.mean(cellCoords, axis=1)
    #target = (-75.0e+3, -50e+3)
    #dist = (cellCenters[:,0]-target[0])**2 + (cellCenters[:,1]-target[1])**2
    #icell = numpy.argmin(dist)

    from pyre.units.time import year
    timePt /= year.value
    timeCell /= year.value
    
    return {'time_pt': timePt,
            'disp': dispPt,
            'vel': velPt, 
            'time_cell': timeCell,
            'cell_centers': cellCenters,
            'stress': stress,
            'cauchy_stress': stressCauchy}

# ======================================================================
prob = "extension"

# Gravity
sim = "%s_gravity" % prob
dataF = getData(sim)

# No-gravity
sim = "%s" % prob
dataI = getData(sim)

import pylab
pylab.subplot(2,2,1)
pylab.plot(dataF['time_pt'], dataF['disp'][:,1], 'r-', 
           dataI['time_pt'], dataI['disp'][:,1], 'b--')

pylab.subplot(2,2,2)
pylab.plot(dataF['time_pt'], dataF['vel'][:,1], 'r-', 
           dataI['time_pt'], dataI['vel'][:,1], 'b--')

pylab.subplot(2,2,3)
icell = 0
pylab.plot(dataF['time_cell'], dataF['stress'][:,icell,0], 'r-', 
           dataF['time_cell'], dataF['cauchy_stress'][:,icell,0], 'g-', 
           dataI['time_cell'], dataI['stress'][:,icell,0], 'b--',
           dataI['time_cell'], dataI['cauchy_stress'][:,icell,0], 'm--')

pylab.subplot(2,2,4)
pylab.plot(dataF['time_cell'], dataF['stress'][:,icell,1], 'r-', 
           dataF['time_cell'], dataF['cauchy_stress'][:,icell,1], 'g-', 
           dataI['time_cell'], dataI['stress'][:,icell,1], 'b--',
           dataI['time_cell'], dataI['cauchy_stress'][:,icell,1], 'm--')

pylab.show()

