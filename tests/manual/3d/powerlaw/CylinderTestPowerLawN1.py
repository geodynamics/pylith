#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/manual/3d/powerlaw/CylinderTestPowerLawN1.py

## @brief Python script to test power-law implementation for steady-state solution
##        of a pressurized cylinder.

import math
import numpy
import h5py
import netCDF4
from pylith.meshio.Xdmf import Xdmf

from cylinderpres_powerlaw_n1_soln import AnalyticalSoln

# ----------------------------------------------------------------------
# Filenames.
h5Prefix = 'output/cylinder_pres_powerlaw_n1_norefstate_'
bqSuffixes = ['b1_q1_', 'b2_q2_']
outputPrefix = 'cylinder_pres_powerlaw_n1_compare_'
cellTypes = ['hex', 'tet']
dispSuffix = '-domain.h5'
stressSuffix = '-viscomat.h5'

# Cylinder parameters.
a = 2000.0
b = 20000.0
h = 2000.0
volCyl = math.pi*h*(b*b - a*a)*0.25

# Tolerances.
volRelTol = 5.0e-3

# ----------------------------------------------------------------------
def stressInvar(stressVec):
    """Compute stress invariants.
    """
    (nsteps, npts, ncomps) = stressVec.shape
    meanStress = ((stressVec[:,:, 0] + stressVec[:,:, 1] + stressVec[:,:, 2])/3.0).reshape(nsteps, npts, 1)
    meanStressVec = meanStress*numpy.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=numpy.float64).reshape(1, 6)
    devStress = stressVec - meanStressVec
    stressInvar = devStress[:,:,0]*devStress[:,:,0] + devStress[:,:,1]*devStress[:,:,1] + devStress[:,:,2]*devStress[:,:,2] \
        + 2.0*(devStress[:,:,3]*devStress[:,:,3] + devStress[:,:,4]*devStress[:,:,4] + devStress[:,:,4]*devStress[:,:,4])
    sqrtInvar = numpy.sqrt(0.5*stressInvar)

    return (sqrtInvar, meanStress)


def cartToCylStress(locs, stressVec):
    """Convert Cartesian stresses to cylindrical.
    """
    angs = numpy.arctan2(locs[:,1], locs[:,0])
    ca = numpy.cos(angs)
    sa = numpy.sin(angs)
    c2a = numpy.cos(2.0*angs)
    s2a = numpy.sin(2.0*angs)
    srr = ca*ca*stressVec[:,:,0] + sa*sa*stressVec[:,:,1] + s2a*stressVec[:,:,3]
    stt = sa*sa*stressVec[:,:,0] + ca*ca*stressVec[:,:,1] - s2a*stressVec[:,:,3]
    srt = -sa*ca*stressVec[:,:,0] + sa*ca*stressVec[:,:,1] + c2a*stressVec[:,:,3]
    stressCyl = numpy.zeros_like(stressVec)
    stressCyl[:,:,0] = srr
    stressCyl[:,:,1] = stt
    stressCyl[:,:,2] = stressVec[:,:,2]
    stressCyl[:,:,3] = srt
    stressCyl[:,:,4] = stressVec[:,:,4]
    stressCyl[:,:,5] = stressVec[:,:,5]

    return stressCyl


def tetVol(vertices, cells):
    """Compute volume of a set of tets.
    """
    (ncells, nvertsPerCell) = cells.shape
    volume = numpy.zeros(ncells, dtype=numpy.float64)
    
    if (nvertsPerCell != 4):
        msg = "Incorrect number of vertices per cell for a tet:  %d" % nvertsPerCell
        raise ValueError(msg)
    cellVerts = vertices[cells, :]
    topRow = numpy.ones(4, dtype=numpy.float64)
    for cellNum in range(ncells):
        verts = cellVerts[cellNum,:]
        vertArr = numpy.vstack((topRow, numpy.transpose(verts)))
        volume[cellNum] = math.fabs(numpy.linalg.det(vertArr)/6.0)

    return volume


def hexVol(vertices, cells):
    """Compute volume of a set of hexes by dividing them into sets of 5 tets.
    """
    (ncells, nvertsPerCell) = cells.shape
    volume = numpy.zeros(ncells, dtype=numpy.float64)
    
    if (nvertsPerCell != 8):
        msg = "Incorrect number of vertices per cell for a hex:  %d" % nvertsPerCell
        raise ValueError(msg)

    c1 = cells[:, [0, 1, 2, 5]]
    c2 = cells[:, [0, 2, 3, 7]]
    c3 = cells[:, [0, 5, 7, 4]]
    c4 = cells[:, [2, 7, 5, 6]]
    c5 = cells[:, [0, 2, 7, 5]]

    v1 = tetVol(vertices, c1)
    v2 = tetVol(vertices, c2)
    v3 = tetVol(vertices, c3)
    v4 = tetVol(vertices, c4)
    v5 = tetVol(vertices, c5)

    volume = v1 + v2 + v3 + v4 + v5

    return volume
        
    
# ----------------------------------------------------------------------
# Loop over cell types.
for cellType in cellTypes:
    for bqSuff in bqSuffixes:
        stressFile = h5Prefix + bqSuff + cellType + stressSuffix
        dispFile = h5Prefix + bqSuff + cellType + dispSuffix
        outputFile = outputPrefix + bqSuff + cellType + '.h5'

        # Read stress info from HDF5 file.
        h5Stress = h5py.File(stressFile, 'r')
        verts = h5Stress['geometry/vertices'][:]
        nverts = verts.shape[0]
        connect = numpy.array(h5Stress['topology/cells'][:], dtype=numpy.int64)
        ncells = connect.shape[0]
        cellCoords = verts[connect, :]
        cellCenters = numpy.mean(cellCoords, axis=1)
        if (cellType == 'tet'):
            cellVols = tetVol(verts, connect)
        else:
            cellVols = hexVol(verts, connect)
        cellVolMesh = numpy.sum(cellVols)
        if (numpy.abs(cellVolMesh - volCyl)/volCyl > volRelTol):
            msg = 'Sum of cell volumes does not match cylinder volume.'
            raise ValueError(msg)
    
        stressNum = h5Stress['cell_fields/cauchy_stress'][:]
        stressCylNum = cartToCylStress(cellCenters, stressNum)
        nsteps = stressNum.shape[0]
        h5Stress.close()
        (stressInvarNum, meanStressNum) = stressInvar(stressNum)

        # Read displacement info from HDF5 file.
        h5Disp = h5py.File(dispFile, 'r')
        time = h5Disp['time'][:,0,0]
        disp = h5Disp['vertex_fields/displacement'][:]
        dt = numpy.diff(time, prepend=0.0)
        d0 = numpy.zeros_like(disp[0,:,:]).reshape(1, nverts, 3)
        dispIncrNum = numpy.diff(disp, axis=0, prepend=d0)
        velocityNum = dispIncrNum/dt.reshape(nsteps, 1, 1)
        h5Disp.close()

        # Compute analytical solution.
        soln = AnalyticalSoln()
        (stressAnl, stressCylAnl) = soln.stress(cellCenters)
        velocityAnl = soln.velocity(verts)
        velocityAnlTiled = numpy.repeat(velocityAnl.reshape(1, nverts, 3), nsteps, axis=0)
        (stressInvarAnl, meanStressAnl) = stressInvar(stressAnl.reshape(1, ncells, 6))
        stressInvarAnlTiled = numpy.repeat(stressInvarAnl.reshape(1, ncells, 1), nsteps, axis=0)
        stressAnlTiled = numpy.repeat(stressAnl.reshape(1, ncells, 6), nsteps, axis=0)
        stressCylAnlTiled = numpy.repeat(stressCylAnl.reshape(1, ncells, 6), nsteps, axis=0)
        meanStressAnlTiled = numpy.repeat(meanStressAnl.reshape(1, ncells, 1), nsteps, axis=0)

        # Compute difference.
        velDiff = velocityAnlTiled - velocityNum
        stressDiff = stressAnlTiled - stressNum
        stressCylDiff = stressCylAnlTiled - stressCylNum
        stressInvarDiff = stressInvarAnlTiled - stressInvarNum.reshape(nsteps, ncells, 1)
        meanStressDiff = meanStressAnlTiled - meanStressNum.reshape(nsteps, ncells, 1)
        stressInvarWtDiff = (cellVols.reshape(ncells, 1, 1)*stressInvarDiff.swapaxes(0,1)/cellVolMesh).swapaxes(0,1)
        meanStressWtDiff = (cellVols.reshape(ncells, 1, 1)*meanStressDiff.swapaxes(0,1)/cellVolMesh).swapaxes(0,1)
        stressInvarDiffSum = numpy.sum(stressInvarDiff)
        meanStressDiffSum = numpy.sum(meanStressDiff)
        stressInvarWtDiffSum = numpy.sum(stressInvarWtDiff)
        meanStressWtDiffSum = numpy.sum(meanStressWtDiff)
        print("Stress invariant summed difference:  %g" % stressInvarDiffSum)
        print("Stress invariant volume-weighted summed difference:  %g" % stressInvarWtDiffSum)
        print("Mean stress summed difference:  %g" % meanStressDiffSum)
        print("Mean stress volume-weighted summed difference:  %g" % meanStressWtDiffSum)

        # Write HDF5 file with results.
        cellDim = 3
        h5Out = h5py.File(outputFile, 'w')
        vertices = h5Out.create_dataset('geometry/vertices', data=verts)
        times = h5Out.create_dataset('time', data=time.reshape(nsteps, 1, 1), maxshape=(None, 1, 1))
        cells = h5Out.create_dataset('topology/cells', data=connect, dtype='d')
        cells.attrs['cell_dim'] = numpy.int32(cellDim)
        velocityAnlOut = h5Out.create_dataset('vertex_fields/velocity_anl',
                                              data=velocityAnlTiled.reshape(nsteps, nverts, 3), maxshape=(None, nverts, 3))
        velocityAnlOut.attrs['vector_field_type'] = 'vector'.encode('ascii')
        velocityNumOut = h5Out.create_dataset('vertex_fields/velocity_num',
                                              data=velocityNum.reshape(nsteps, nverts, 3), maxshape=(None, nverts, 3))
        velocityNumOut.attrs['vector_field_type'] = 'vector'.encode('ascii')
        velocityAnlMinusNum = h5Out.create_dataset('vertex_fields/velocity_anl_minus_num',
                                                   data=velDiff.reshape(nsteps, nverts, 3), maxshape=(None, nverts, 3))
        velocityAnlMinusNum.attrs['vector_field_type'] = 'vector'.encode('ascii')
        stressVecAnl = h5Out.create_dataset('cell_fields/stress_anl',
                                            data=stressAnlTiled.reshape(nsteps, ncells, 6), maxshape=(None, ncells, 6))
        stressVecAnl.attrs['vector_field_type'] = 'tensor'.encode('ascii')
        stressVecCylAnl = h5Out.create_dataset('cell_fields/stress_cylindrical_anl',
                                               data=stressCylAnlTiled.reshape(nsteps, ncells, 6), maxshape=(None, ncells, 6))
        stressVecCylAnl.attrs['vector_field_type'] = 'tensor'.encode('ascii')
        stressVecNum = h5Out.create_dataset('cell_fields/stress_num',
                                            data=stressNum.reshape(nsteps, ncells, 6), maxshape=(None, ncells, 6))
        stressVecNum.attrs['vector_field_type'] = 'tensor'.encode('ascii')
        stressVecCylNum = h5Out.create_dataset('cell_fields/stress_cylindrical_num',
                                               data=stressCylNum.reshape(nsteps, ncells, 6), maxshape=(None, ncells, 6))
        stressVecCylNum.attrs['vector_field_type'] = 'tensor'.encode('ascii')
        stressAnlMinusNum = h5Out.create_dataset('cell_fields/stress_anl_minus_num',
                                                 data=stressDiff.reshape(nsteps, ncells, 6), maxshape=(None, ncells, 6))
        stressAnlMinusNum.attrs['vector_field_type'] = 'tensor'.encode('ascii')
        stressAnlMinusNumCyl = h5Out.create_dataset('cell_fields/stress_anl_minus_num_cyl',
                                                    data=stressCylDiff.reshape(nsteps, ncells, 6), maxshape=(None, ncells, 6))
        stressAnlMinusNumCyl.attrs['vector_field_type'] = 'tensor'.encode('ascii')
        stressInvar2Anl = h5Out.create_dataset('cell_fields/stress_invar_anl',
                                               data=stressInvarAnlTiled.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
        stressInvar2Anl.attrs['vector_field_type'] = 'scalar'.encode('ascii')
        stressInvar2Num = h5Out.create_dataset('cell_fields/stress_invar_num',
                                               data=stressInvarNum.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
        stressInvar2Num.attrs['vector_field_type'] = 'scalar'.encode('ascii')
        stressInvar2Diff = h5Out.create_dataset('cell_fields/stress_invar_anl_minus_num',
                                                data=stressInvarDiff.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
        stressInvar2Diff.attrs['vector_field_type'] = 'scalar'.encode('ascii')
        stressInvar2WtDiff = h5Out.create_dataset('cell_fields/stress_invar_anl_minus_num_vol_weighted',
                                                  data=stressInvarWtDiff.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
        stressInvar2WtDiff.attrs['vector_field_type'] = 'scalar'.encode('ascii')
        meanStress2Anl = h5Out.create_dataset('cell_fields/mean_stress_anl',
                                              data=meanStressAnlTiled.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
        meanStress2Anl.attrs['vector_field_type'] = 'scalar'.encode('ascii')
        meanStress2Num = h5Out.create_dataset('cell_fields/mean_stress_num',
                                              data=meanStressNum.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
        meanStress2Num.attrs['vector_field_type'] = 'scalar'.encode('ascii')
        meanStress2Diff = h5Out.create_dataset('cell_fields/mean_stress_anl_minus_num',
                                               data=meanStressDiff.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
        meanStress2Diff.attrs['vector_field_type'] = 'scalar'.encode('ascii')
        meanStress2WtDiff = h5Out.create_dataset('cell_fields/mean_stress_anl_minus_num_vol_weighted',
                                                 data=meanStressWtDiff.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
        meanStress2WtDiff.attrs['vector_field_type'] = 'scalar'.encode('ascii')
        cellVolRepeat = numpy.repeat(cellVols.reshape(1, ncells, 1), nsteps, axis=0)
        cellVol = h5Out.create_dataset('cell_fields/cell_volume',
                                       data=cellVolRepeat.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
        cellVol.attrs['vector_field_type'] = 'scalar'.encode('ascii')
        h5Out.close()

        xdmfWriter = Xdmf()
        xdmfWriter.write(outputFile)
    

# End of file 
