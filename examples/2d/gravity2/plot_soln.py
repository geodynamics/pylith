import h5py
import numpy

def getData(sim):

    filename = "output/%s-groundsurf.h5" % sim
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

    filename = "output/%s-mantle.h5" % sim
    h5 = h5py.File(filename, "r")
    vertices = h5['geometry/vertices'][:]
    cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
    vstrain = h5['cell_fields/viscous_strain'][:]
    stress = h5['cell_fields/cauchy_stress'][-1,:,:]
    timeCell = h5['time'][:,0,0]
    h5.close()

    cellCoords = vertices[cells,:]
    cellCenters = numpy.mean(cellCoords, axis=1)
    target = (-75.0e+3, -50e+3)
    dist = (cellCenters[:,0]-target[0])**2 + (cellCenters[:,1]-target[1])**2
    icell = numpy.argmin(dist)
    vstrainCell = vstrain[:,icell,:]

    from pyre.units.time import year
    timePt /= year.value
    timeCell /= year.value
    
    return {'time_pt': timePt,
            'disp': dispPt,
            'vel': velPt, 
            'time_cell': timeCell,
            'vstrain': vstrainCell, 
            'cell_centers': cellCenters,
            'stress': stress}

# ======================================================================
geom = "topo"

# Finite strain
sim = "%s_inf" % geom
dataF = getData(sim)

# Infinitesimal strain
sim = "%s_inf" % geom
dataI = getData(sim)

import pylab
pylab.subplot(2,2,1)
pylab.plot(dataF['time_pt'], dataF['disp'][:,1], 'r-', 
           dataI['time_pt'], dataI['disp'][:,1], 'b--')

pylab.subplot(2,2,2)
pylab.plot(dataF['time_pt'], dataF['vel'][:,1], 'r-', 
           dataI['time_pt'], dataI['vel'][:,1], 'b--')

pylab.subplot(2,2,3)
pylab.plot(dataF['time_cell'], dataF['vstrain'][:,:], 'r-', 
           dataI['time_cell'], dataI['vstrain'][:,:], 'b--')

pylab.subplot(2,2,4)
def stressDiff(data):
    return (data['stress'][:,1] - data['stress'][:,0]) / data['stress'][:,1]
stressDiffF = stressDiff(dataF)
stressDiffI = stressDiff(dataI)
pylab.plot(dataF['stress'][:,0], dataF['cell_centers'][:,1], 'r.', 
           dataI['stress'][:,0], dataI['cell_centers'][:,1], 'bx')

pylab.show()

