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

    target = (0.0, 0.0)
    dist = (vertices[:,0]-target[0])**2 + (vertices[:,1]-target[1])**2
    ipt = numpy.argmin(dist)
    dispPt = disp[:,ipt,:]
    velPt = vel[:,ipt,:]

    filename = "output/%s-mantle.h5" % sim
    h5 = h5py.File(filename, "r")
    vertices = h5['geometry/vertices'][:]
    cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
    vstrain = h5['cell_fields/viscous_strain'][:]
    stress = h5['cell_fields/cauchy_stress'][:,:,:]
    timeCell = h5['time'][:,0,0]
    h5.close()

    cellCoords = vertices[cells,:]
    cellCenters = numpy.mean(cellCoords, axis=1)
    target = (-50.0e+3, -50e+3)
    dist = (cellCenters[:,0]-target[0])**2 + (cellCenters[:,1]-target[1])**2
    icell = numpy.argmin(dist)
    vstrainCell = vstrain[:,icell,:]
    stressCell = stress[:,icell,:]

    from pyre.units.time import year
    timePt /= year.value
    timeCell /= year.value
    
    return {'time_pt': timePt,
            'disp': dispPt,
            'vel': velPt, 
            'time_cell': timeCell,
            'vstrain': vstrainCell, 
            'cell_centers': cellCenters,
            'stress': stressCell}

# ======================================================================
sim = "postseismic_vardensity"
dataA = getData(sim)

sim = "postseismic_finstrain"
dataB = getData(sim)

import pylab
pylab.subplot(2,2,1)
pylab.plot(dataA['time_pt'], dataA['disp'][:,1], 'r-',
           dataB['time_pt'], dataB['disp'][:,1], 'b--')
pylab.xlabel('Time (year)')
pylab.ylabel('Displacement (m)')

pylab.subplot(2,2,2)
pylab.plot(dataA['time_pt'], dataA['vel'][:,:], 'r-',
           dataB['time_pt'], dataB['vel'][:,:], 'b--')
pylab.xlabel('Time (year)')
pylab.ylabel('Velocity (m/s)')

pylab.subplot(2,2,3)
pylab.plot(dataA['time_cell'], dataA['vstrain'][:,:], 'r-',
           dataB['time_cell'], dataB['vstrain'][:,:], 'b--',)
pylab.xlabel('Time (year)')
pylab.ylabel('Viscous Strain')

pylab.subplot(2,2,4)
pylab.plot(dataA['time_cell'], dataA['stress'][:,2], 'r-',
           dataB['time_cell'], dataB['stress'][:,2], 'b--',)
pylab.xlabel('Time (year)')
pylab.ylabel('Stress_xy')

pylab.show()

