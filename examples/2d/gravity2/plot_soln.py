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
sim = "initialstress"
data = getData(sim)

import pylab
pylab.subplot(2,2,1)
pylab.plot(data['time_pt'], data['disp'][:,1], 'r-')
pylab.xlabel('Time (year)')
pylab.ylabel('Displacement (m)')

pylab.subplot(2,2,2)
pylab.plot(data['time_pt'], data['vel'][:,1], 'r-')
pylab.xlabel('Time (year)')
pylab.ylabel('Velocity (m/s)')

pylab.subplot(2,2,3)
pylab.plot(data['time_cell'], data['vstrain'][:,:], 'r-')
pylab.xlabel('Time (year)')
pylab.ylabel('Viscous Strain')

pylab.show()

