import h5py
import numpy

def getData(sim):

    filename = "%s-domain.h5" % sim
    h5 = h5py.File(filename, "r")
    vertices = h5['geometry/vertices'][:]
    disp = h5['vertex_fields/displacement'][0,:,:]
    h5.close()

    filename = "%s-elastic.h5" % sim
    h5 = h5py.File(filename, "r")
    vertices = h5['geometry/vertices'][:]
    cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)
    strain = h5['cell_fields/total_strain'][0,:,:]
    stress = h5['cell_fields/stress'][0,:,:]
    stressCauchy = h5['cell_fields/cauchy_stress'][0,:,:]
    h5.close()

    cellCoords = vertices[cells,:]
    cellCenters = numpy.mean(cellCoords, axis=1)
    #target = (-75.0e+3, -50e+3)
    #dist = (cellCenters[:,0]-target[0])**2 + (cellCenters[:,1]-target[1])**2
    #icell = numpy.argmin(dist)

    return {'disp': disp,
            'cell_centers': cellCenters,
            'strain': strain,
            'stress': stress,
            'cauchy_stress': stressCauchy}

# ======================================================================
prob = "compression_gravity"

# Finite strain
sim = "%s" % prob
dataF = getData(sim)

# Infinitesimal
sim = "%s_inf" % prob
dataI = getData(sim)

def strainEnergy(data):
    return 0.5*0.5*numpy.sum(data['stress']*data['strain'])*2.0e+3**2

def work(data, finite=False):
    g = 9.80665
    density = 2500.0
    disp = numpy.abs(data['disp'][1,1])
    if finite:
        l = 2.0e+3
        from math import log
        disp = -l*log(1.0-disp/l)
    return 2*g*density*0.5*disp*2.0e+3**2

print "Infinitesimal"
print "Disp",dataI['disp']
print "Strain",dataI['strain']
print "Stress",dataI['stress']
print "Strain energy",strainEnergy(dataI)
print "Work",work(dataI)

print "Finite"
print "Disp",dataF['disp']
print "Strain",dataF['strain']
print "Stress",dataF['stress']
print "Cauchy Stress",dataF['cauchy_stress']
print "Strain energy",strainEnergy(dataF)
print "Work",work(dataF, finite=True)
