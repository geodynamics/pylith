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

## @file tests/manual/2d/powerlaw/CylinderTestElastic.py

## @brief Python script to test elastic implementation for a pressurized cylinder.

import math
import numpy
import h5py
import netCDF4
from pylith.meshio.Xdmf import Xdmf

from cylinderpres_elastic_soln import AnalyticalSoln

# ----------------------------------------------------------------------
# Filenames.
h5Prefix = 'output/cylinder_pres_elastic_norefstate_'
bqSuffixes = ['b1_q1', 'b2_q2', 'b3_q3', 'b4_q4']
meshPrefix = 'meshes/mesh_cylinder'
meshResolutions = ['_', '_fine_']
outputPrefix = 'cylinder_pres_elastic_compare'
cellTypes = ['quad', 'tri']
dispSuffix = '-domain.h5'
stressSuffix = '-elastic.h5'

# Cylinder parameters.
a = 2000.0
b = 20000.0
areaCyl = math.pi*(b*b - a*a)*0.25

areaRelTol = 5.0e-3

# ----------------------------------------------------------------------
def stressInvar(stressVec):
    """Compute stress invariants.
    """
    (nsteps, npts, ncomps) = stressVec.shape
    meanStress = ((stressVec[:,:, 0] + stressVec[:,:, 1] + stressVec[:,:, 2])/3.0).reshape(nsteps, npts, 1)
    meanStressVec = meanStress*numpy.array([1.0, 1.0, 1.0, 0.0], dtype=numpy.float64).reshape(1, 4)
    devStress = stressVec - meanStressVec
    stressInvar = devStress[:,:,0]*devStress[:,:,0] + devStress[:,:,1]*devStress[:,:,1] + devStress[:,:,2]*devStress[:,:,2] \
        + 2.0*devStress[:,:,3]*devStress[:,:,3]
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

    return stressCyl


def triArea(vertices, cells):
    """Compute areas of a set of tris.
    """
    (ncells, nvertsPerCell) = cells.shape
    area = numpy.zeros(ncells, dtype=numpy.float64)
    
    if (nvertsPerCell != 3):
        msg = "Incorrect number of vertices per cell for a tri:  %d" % nvertsPerCell
        raise ValueError(msg)
    cellVerts = vertices[cells, :]
    topRow = numpy.ones(3, dtype=numpy.float64)
    for cellNum in range(ncells):
        verts = cellVerts[cellNum,:]
        vertArr = numpy.vstack((topRow, numpy.transpose(verts)))
        area[cellNum] = math.fabs(numpy.linalg.det(vertArr)/2.0)

    return area


def quadArea(vertices, cells):
    """Compute areas of a set of quads by dividing them into sets of 2 tris.
    """
    (ncells, nvertsPerCell) = cells.shape
    area = numpy.zeros(ncells, dtype=numpy.float64)
    
    if (nvertsPerCell != 4):
        msg = "Incorrect number of vertices per cell for a quad:  %d" % nvertsPerCell
        raise ValueError(msg)

    c1 = cells[:, [0, 1, 2]]
    c2 = cells[:, [0, 2, 3]]

    s1 = triArea(vertices, c1)
    s2 = triArea(vertices, c2)

    area = s1 + s2

    return area
        
    
# ----------------------------------------------------------------------
# Loop over cell types.
for cellType in cellTypes:
    for meshResolution in meshResolutions:
        for bqSuff in bqSuffixes:
            rootName = bqSuff + meshResolution + cellType
            stressFile = h5Prefix + rootName + stressSuffix
            dispFile = h5Prefix + rootName + dispSuffix
            outputFile = outputPrefix + rootName + '.h5'

            # Read stress info from HDF5 file.
            h5Stress = h5py.File(stressFile, 'r')
            verts = h5Stress['geometry/vertices'][:]
            nverts = verts.shape[0]
            connect = numpy.array(h5Stress['topology/cells'][:], dtype=numpy.int64)
            ncells = connect.shape[0]
            cellCoords = verts[connect, :]
            cellCenters = numpy.mean(cellCoords, axis=1)
            if (cellType == 'tri'):
                cellAreas = triArea(verts, connect)
            else:
                cellAreas = quadArea(verts, connect)
            cellAreaMesh = numpy.sum(cellAreas)
            if (numpy.abs(cellAreaMesh - areaCyl)/areaCyl > areaRelTol):
                msg = 'Sum of cell areas does not match cylinder area.'
                raise ValueError(msg)
    
            stressNum = h5Stress['cell_fields/cauchy_stress'][:]
            stressCylNum = cartToCylStress(cellCenters, stressNum)
            nsteps = stressNum.shape[0]
            h5Stress.close()
            (stressInvarNum, meanStressNum) = stressInvar(stressNum)

            # Read displacement info from HDF5 file.
            h5Disp = h5py.File(dispFile, 'r')
            time = h5Disp['time'][:,0,0]
            dispNum = h5Disp['vertex_fields/displacement'][:]
            h5Disp.close()

            # Compute analytical solution.
            soln = AnalyticalSoln()
            (stressAnl, stressCylAnl) = soln.stress(cellCenters)
            dispAnl = soln.displacement(verts)
            dispAnlTiled = numpy.repeat(dispAnl.reshape(1, nverts, 2), nsteps, axis=0)
            (stressInvarAnl, meanStressAnl) = stressInvar(stressAnl.reshape(1, ncells, 4))
            stressInvarAnlTiled = numpy.repeat(stressInvarAnl.reshape(1, ncells, 1), nsteps, axis=0)
            stressAnlTiled = numpy.repeat(stressAnl.reshape(1, ncells, 4), nsteps, axis=0)
            stressCylAnlTiled = numpy.repeat(stressCylAnl.reshape(1, ncells, 4), nsteps, axis=0)
            meanStressAnlTiled = numpy.repeat(meanStressAnl.reshape(1, ncells, 1), nsteps, axis=0)

            # Compute difference.
            dispDiff = dispAnlTiled - dispNum
            stressDiff = stressAnlTiled - stressNum
            stressCylDiff = stressCylAnlTiled - stressCylNum
            stressInvarDiff = stressInvarAnlTiled - stressInvarNum.reshape(nsteps, ncells, 1)
            meanStressDiff = meanStressAnlTiled - meanStressNum.reshape(nsteps, ncells, 1)
            stressInvarWtDiff = (cellAreas.reshape(ncells, 1, 1)*stressInvarDiff.swapaxes(0,1)/cellAreaMesh).swapaxes(0,1)
            meanStressWtDiff = (cellAreas.reshape(ncells, 1, 1)*meanStressDiff.swapaxes(0,1)/cellAreaMesh).swapaxes(0,1)
            stressInvarDiffSum = numpy.sum(numpy.sqrt(stressInvarDiff*stressInvarDiff))
            meanStressDiffSum = numpy.sum(numpy.sqrt(meanStressDiff*meanStressDiff))
            stressInvarWtDiffSum = numpy.sum(numpy.sqrt(stressInvarWtDiff*stressInvarWtDiff))
            meanStressWtDiffSum = numpy.sum(numpy.sqrt(meanStressWtDiff*meanStressWtDiff))
            print("Model:  %s" % rootName)
            print("Stress invariant summed difference:  %g" % stressInvarDiffSum)
            print("Stress invariant volume-weighted summed difference:  %g" % stressInvarWtDiffSum)
            print("Mean stress summed difference:  %g" % meanStressDiffSum)
            print("Mean stress volume-weighted summed difference:  %g" % meanStressWtDiffSum)
            print("")

            # Write HDF5 file with results.
            cellDim = 2
            h5Out = h5py.File(outputFile, 'w')
            vertices = h5Out.create_dataset('geometry/vertices', data=verts)
            times = h5Out.create_dataset('time', data=time.reshape(nsteps, 1, 1), maxshape=(None, 1, 1))
            cells = h5Out.create_dataset('topology/cells', data=connect, dtype='d')
            cells.attrs['cell_dim'] = numpy.int32(cellDim)
            dispAnlOut = h5Out.create_dataset('vertex_fields/displacement_anl',
                                              data=dispAnlTiled.reshape(nsteps, nverts, 2), maxshape=(None, nverts, 2))
            dispAnlOut.attrs['vector_field_type'] = 'vector'.encode('ascii')
            dispNumOut = h5Out.create_dataset('vertex_fields/displacement_num',
                                              data=dispNum.reshape(nsteps, nverts, 2), maxshape=(None, nverts, 2))
            dispNumOut.attrs['vector_field_type'] = 'vector'.encode('ascii')
            dispAnlMinusNum = h5Out.create_dataset('vertex_fields/displacement_anl_minus_num',
                                                   data=dispDiff.reshape(nsteps, nverts, 2), maxshape=(None, nverts, 2))
            dispAnlMinusNum.attrs['vector_field_type'] = 'vector'.encode('ascii')
            stressVecAnl = h5Out.create_dataset('cell_fields/stress_anl',
                                                data=stressAnlTiled.reshape(nsteps, ncells, 4), maxshape=(None, ncells, 4))
            stressVecAnl.attrs['vector_field_type'] = 'other'.encode('ascii')
            stressVecCylAnl = h5Out.create_dataset('cell_fields/stress_cylindrical_anl',
                                                   data=stressCylAnlTiled.reshape(nsteps, ncells, 4), maxshape=(None, ncells, 4))
            stressVecCylAnl.attrs['vector_field_type'] = 'other'.encode('ascii')
            stressVecNum = h5Out.create_dataset('cell_fields/stress_num',
                                                data=stressNum.reshape(nsteps, ncells, 4), maxshape=(None, ncells, 4))
            stressVecNum.attrs['vector_field_type'] = 'other'.encode('ascii')
            stressVecCylNum = h5Out.create_dataset('cell_fields/stress_cylindrical_num',
                                                   data=stressCylNum.reshape(nsteps, ncells, 4), maxshape=(None, ncells, 4))
            stressVecCylNum.attrs['vector_field_type'] = 'other'.encode('ascii')
            stressAnlMinusNum = h5Out.create_dataset('cell_fields/stress_anl_minus_num',
                                                     data=stressDiff.reshape(nsteps, ncells, 4), maxshape=(None, ncells, 4))
            stressAnlMinusNum.attrs['vector_field_type'] = 'other'.encode('ascii')
            stressAnlMinusNumCyl = h5Out.create_dataset('cell_fields/stress_anl_minus_num_cyl',
                                                        data=stressCylDiff.reshape(nsteps, ncells, 4), maxshape=(None, ncells, 4))
            stressAnlMinusNumCyl.attrs['vector_field_type'] = 'other'.encode('ascii')
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
            cellAreaRepeat = numpy.repeat(cellAreas.reshape(1, ncells, 1), nsteps, axis=0)
            cellArea = h5Out.create_dataset('cell_fields/cell_area',
                                            data=cellAreaRepeat.reshape(nsteps, ncells, 1), maxshape=(None, ncells, 1))
            cellArea.attrs['vector_field_type'] = 'scalar'.encode('ascii')
            h5Out.close()

            xdmfWriter = Xdmf()
            xdmfWriter.write(outputFile)
    

# End of file 
