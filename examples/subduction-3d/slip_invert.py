#!/usr/bin/env nemesis
# -*- Python -*- (syntax highlighting)
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# Simple inversion script that uses PyLith-generated Green's functions.
#
# This script is used as part of example step07. To perform the inversion you
# must have already:
# 1.  Run example step06.
# 2.  Run the make_synthetic_gpsdisp.py script to generate synthetic data.
# 3.  Generated the step07 Green's functions (step07a and step07b).
#
# Once you have performed the steps above, you can run this script. The
# parameters are defined in slip_invert.cfg.
# Run this script as follows:
# ./slip_invert.py
#

import math
import numpy
import sys
import os
from pythia.pyre.units.angle import degree
import h5py

from pythia.pyre.applications.Script import Script as Application


class SlipInvert(Application):
    """Python application to perform a linear inversion for slip using
    PyLith-generated Green's functions.
    """

    # \b Properties
    # @li \b data_file File with displacements, locations, and stdDev.
    # @li \b rake Assumed rake direction.
    # @li \b gfimpulses_ll_file HDF5 file with left-lateral GF impulses.
    # @li \b gfimpulses_ud_file HDF5 file with updip GF impulses.
    # @li \b gfresponses_ll_file HDF5 file with left-lateral GF responses.
    # @li \b gfresponses_ud_file HDF5 file with updip GF responses.
    # @li \b a_priori_value A priori value for parameters.
    # @li \b penalty_weight_vals List of penalty weights.
    # @li \b data_scale Scaling factor to apply to data and stdDev.
    # @li \b result_summary_file Text file summarizing inversion results.
    # @li \b slip_output_file Output file with slip results.
    # @li \b displacement_output_file Output file with inversion results.
    ##
    # \b Facilities
    # @li None

    import pythia.pyre.inventory

    dataFile = pythia.pyre.inventory.str("data_file", default="data.txt")
    dataFile.meta['tip'] = "File with displ., locations, and stdDev."

    rake = pythia.pyre.inventory.dimensional("rake", default=90.0 * degree)
    rake.meta['tip'] = "Assumed rake angle."

    gfImpulsesLlFile = pythia.pyre.inventory.str("gfimpulses_ll_file", default="gfimpulse_ll.h5")
    gfImpulsesLlFile.meta['tip'] = "HDF5 file with left-lateral GF impulses."

    gfImpulsesUdFile = pythia.pyre.inventory.str("gfimpulses_ud_file", default="gfimpulse_ud.h5")
    gfImpulsesUdFile.meta['tip'] = "HDF5 file with updip GF impulses."

    gfResponsesLlFile = pythia.pyre.inventory.str("gfresponses_ll_file", default="gfresponse_ll.h5")
    gfResponsesLlFile.meta['tip'] = "HDF5 file with left-lateral GF responses."

    gfResponsesUdFile = pythia.pyre.inventory.str("gfresponses_ud_file", default="gfresponse_ud.h5")
    gfResponsesUdFile.meta['tip'] = "HDF5 file with updip GF responses."

    aPrioriValue = pythia.pyre.inventory.float("a_priori_value", default=0.0)
    aPrioriValue.meta['tip'] = "A priori value for parameters."

    penaltyWeightVals = pythia.pyre.inventory.list("penalty_weight_vals", default=[0.1, 0.5, 1.0, 5.0, 10.0])
    penaltyWeightVals.meta['tip'] = "List of penalty weights."

    dataScale = pythia.pyre.inventory.float("data_scale", default=1.0)
    dataScale.meta['tip'] = "Scaling factor to apply to data and stdDev."

    resultSummaryFile = pythia.pyre.inventory.str("result_summary_file", default='result_summary.txt')
    resultSummaryFile.meta['tip'] = "Text file summarizing inversion results."

    slipOutputFile = pythia.pyre.inventory.str("slip_output_file", default='predicted_slip.h5')
    slipOutputFile.meta['tip'] = "HDF5 file with predicted slip results."

    displacementOutputFile = pythia.pyre.inventory.str("displacement_output_file", default='predicted_displacement.h5')
    displacementOutputFile.meta['tip'] = "HDF5 file with predicted displacements."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="slip_invert"):
        Application.__init__(self, name)

        self.dataCoords = None
        self.dataVals = None
        self.dataCov = None
        self.dataNames = []
        self.numDataPoints = 0
        self.numDesignRows = 0

        self.numFaultVerts = 0
        self.numFaultCells = 0
        self.faultVertCoords = None
        self.faultCells = None

        self.numImpulses = 0
        self.impulseInds = None
        self.impulseCoords = None

        self.summaryHead = 'Penalty-weight\tData-residual\t' + \
                           'Weighted-data-residual\tPenalty-residual\t' + \
                           'Weighted-penalty-residual\tTotal-residual\t' + \
                           'Total-weighted-residual'
        self.numSummaryCols = 7

        self.design = None

        return

    def main(self):
        self.readData()
        self.readGreens()
        self.runInversions()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        Application._configure(self)
        self.penaltyWeights = numpy.array(self.penaltyWeightVals, dtype=numpy.float64)
        self.numPenaltyWeights = self.penaltyWeights.shape[0]

        # Left-lateral and updip components from assumed rake.
        self.llComp = math.cos(self.rake.value)
        self.udComp = math.sin(self.rake.value)

        return

    def runInversions(self):
        """Function to run inversions using a range of penalty parameters.
        """
        print("Running inversions:")
        sys.stdout.flush()

        # Open output files.
        d = h5py.File(self.displacementOutputFile, 'w')
        s = h5py.File(self.slipOutputFile, 'w')

        # Write fault mesh and time info.
        summaryInfo = numpy.zeros((self.numPenaltyWeights, self.numSummaryCols),
                                  dtype=numpy.float64)
        cellDimF = 2
        timesF = self.penaltyWeights.reshape(self.numPenaltyWeights, 1, 1)
        vertsF = s.create_dataset('geometry/vertices', data=self.faultVertCoords)
        timesF = s.create_dataset('time', data=timesF, maxshape=(None, 1, 1))
        topoF = s.create_dataset('topology/cells', data=self.faultCells, dtype='d')
        topoF.attrs['cell_dim'] = numpy.int32(cellDimF)
        slipVec = numpy.zeros((self.numPenaltyWeights, self.numFaultVerts, 3),
                              dtype=numpy.float64)
        slipAlongRake = numpy.zeros((self.numPenaltyWeights, self.numFaultVerts, 1),
                                    dtype=numpy.float64)

        # Write data mesh and time info.
        cellDimD = 0
        topolD = numpy.arange(self.numDataPoints, dtype=numpy.int64).reshape(self.numDataPoints, 1)
        timesD = self.penaltyWeights.reshape(self.numPenaltyWeights, 1, 1)
        vertsD = d.create_dataset('geometry/vertices', data=self.dataCoords)
        timesD = d.create_dataset('time', data=timesD, maxshape=(None, 1, 1))
        topoD = d.create_dataset('topology/cells', data=topolD, dtype='d')
        topoD.attrs['cell_dim'] = numpy.int32(cellDimD)

        predictedDisp = numpy.zeros((self.numPenaltyWeights, self.numDataPoints, 3), dtype=numpy.float64)

        # Rescale equations using data standard deviations.
        dataStdDev = numpy.sqrt(self.dataCov)
        dataStdDevInvDiag = numpy.diag(1.0 / dataStdDev)
        dataScaledDesign = numpy.dot(dataStdDevInvDiag, self.design)
        dataScaledVals = numpy.dot(dataStdDevInvDiag, self.dataVals)

        # Create a priori parameter vector.
        paramVec = self.aPrioriValue * numpy.ones(self.numImpulses, dtype=numpy.float64)

        summFmt = '%g' + 6 * '\t%e' + '\n'

        # Regularization array is just the identity matrix.
        regArray = numpy.identity(self.numImpulses, dtype=numpy.float64)

        # Loop over inversions.
        for invNum in range(self.numPenaltyWeights):
            penWeight = self.penaltyWeights[invNum]
            print("  Working on inversion %d, penalty weight = %g" % (invNum, penWeight))
            sys.stdout.flush()
            paramScaledDesign = penWeight * regArray
            paramScaledData = penWeight * paramVec
            designMat = numpy.vstack((dataScaledDesign, paramScaledDesign))
            dataVec = numpy.hstack((dataScaledVals, paramScaledData))
            designMatTrans = numpy.transpose(designMat)
            genInv = numpy.dot(numpy.linalg.inv(numpy.dot(designMatTrans, designMat)), designMatTrans)
            solution = numpy.dot(genInv, dataVec)

            # Compute residuals, etc.
            predicted = numpy.dot(self.design, solution)
            dataResidual = self.dataVals - predicted
            dataWeightResidual = numpy.dot(dataStdDevInvDiag, dataResidual)
            dataResidualNorm = numpy.linalg.norm(dataResidual)
            dataWeightResidualNorm = numpy.linalg.norm(dataWeightResidual)
            penalty = numpy.dot(regArray, solution)
            penaltyResidual = paramVec - penalty
            penaltyWeightResidual = penWeight * penaltyResidual
            penaltyResidualNorm = numpy.linalg.norm(penaltyResidual)
            penaltyWeightResidualNorm = numpy.linalg.norm(penaltyWeightResidual)
            totalResidualNorm = dataResidualNorm + penaltyResidualNorm
            totalWeightResidualNorm = dataWeightResidualNorm + penaltyWeightResidualNorm

            summaryInfo[invNum, 0] = penWeight
            summaryInfo[invNum, 1] = dataResidualNorm
            summaryInfo[invNum, 2] = dataWeightResidualNorm
            summaryInfo[invNum, 3] = penaltyResidualNorm
            summaryInfo[invNum, 4] = penaltyWeightResidualNorm
            summaryInfo[invNum, 5] = totalResidualNorm
            summaryInfo[invNum, 6] = totalWeightResidualNorm

            slipAlongRake[invNum, self.impulseInds, 0] = solution
            slipVec[invNum, self.impulseInds, 0] = self.llComp * solution
            slipVec[invNum, self.impulseInds, 1] = self.udComp * solution
            predictedDisp[invNum,:,:] = predicted.reshape(self.numDataPoints, 3,
                                                            order='F')

            print("    Data residual:              %e" % dataResidualNorm)
            print("    Weighted data residual:     %e" % dataWeightResidualNorm)
            print("    Penalty residual:           %e" % penaltyResidualNorm)
            print("    Weighted penalty residual:  %e" % penaltyWeightResidualNorm)
            print("    Total residual:             %e" % totalResidualNorm)
            print("    Weighted total residual:    %e" % totalWeightResidualNorm)
            sys.stdout.flush()

        numpy.savetxt(self.resultSummaryFile, summaryInfo, delimiter='\t',
                      header=self.summaryHead)

        # Write results to HDF5 files.
        rakeSlip = s.create_dataset('vertex_fields/rake_slip', data=slipAlongRake)
        rakeSlip.attrs['vector_field_type'] = 'scalar'
        vecSlip = s.create_dataset('vertex_fields/slip_vector', data=slipVec)
        vecSlip.attrs['vector_field_type'] = 'vector'
        vecDisp = d.create_dataset('vertex_fields/disp_vec', data=predictedDisp)
        vecDisp.attrs['vector_field_type'] = 'vector'

        s.close()
        d.close()

        from pylith.meshio.Xdmf import Xdmf
        xdmfWriter = Xdmf()
        xdmfWriter.write(self.slipOutputFile)
        xdmfWriter.write(self.displacementOutputFile)
        return

    def readGreens(self):
        """Function to read impulse and response info from PyLith output files.
        """
        print("Reading Green's functions:")
        sys.stdout.flush()

        # Read impulses.
        print("  Reading left-lateral impulses:")
        sys.stdout.flush()
        impulsesLl = h5py.File(self.gfImpulsesLlFile, 'r')
        self.faultVertCoords = impulsesLl['geometry/vertices'][:]
        self.numFaultVerts = self.faultVertCoords.shape[0]
        self.faultCells = numpy.array(impulsesLl['topology/cells'][:],
                                      dtype=numpy.int)
        self.numFaultCells = self.faultCells.shape[0]
        llSlip = impulsesLl['vertex_fields/slip'][:,:, 0]
        llImpInds = numpy.nonzero(llSlip != 0.0)
        self.impulseCoords = self.faultVertCoords[llImpInds[1]]
        self.numImpulses = self.impulseCoords.shape[0]

        print("    Number of fault vertices:     %d" % self.numFaultVerts)
        print("    Number of impulses:           %d" % self.numImpulses)

        (distances, self.impulseInds) = self.matchCoords(self.faultVertCoords,
                                                         self.impulseCoords)
        impulsesLl.close()

        print("  Reading updip impulses:")
        sys.stdout.flush()
        impulsesUd = h5py.File(self.gfImpulsesUdFile, 'r')
        udCoords = impulsesUd['geometry/vertices'][:]
        udSlip = impulsesUd['vertex_fields/slip'][:,:, 1]
        udImpInds = numpy.nonzero(udSlip != 0.0)
        numUdImpulses = udImpInds[0].shape[0]
        udCoordsUsed = udCoords[udImpInds[1]]
        udSlipUsed = udSlip[udImpInds[0], udImpInds[1]]
        (distances, udCoordInds) = self.matchCoords(udCoordsUsed,
                                                    self.impulseCoords)
        udCoordsUsed = udCoordsUsed[udCoordInds,:]
        impulsesUd.close()

        # Read responses.
        print("  Reading left-lateral responses:")
        sys.stdout.flush()
        responseLl = h5py.File(self.gfResponsesLlFile, 'r')
        llResponseCoords = responseLl['geometry/vertices'][:]
        llResponseVals = responseLl['vertex_fields/displacement'][:]

        (distances, llDataInds) = self.matchCoords(llResponseCoords,
                                                   self.dataCoords)
        llResponsesEast = llResponseVals[:, llDataInds, 0]
        llResponsesNorth = llResponseVals[:, llDataInds, 1]
        llResponsesUp = llResponseVals[:, llDataInds, 2]
        responseLl.close()

        print("  Reading updip responses:")
        sys.stdout.flush()
        responseUd = h5py.File(self.gfResponsesUdFile, 'r')
        udResponseCoords = responseUd['geometry/vertices'][:]
        responseUdVals = responseUd['vertex_fields/displacement'][:]
        udResponseVals = responseUdVals[udCoordInds,:,:]

        (distances, udDataInds) = self.matchCoords(udResponseCoords,
                                                   self.dataCoords)
        udResponsesEast = udResponseVals[:, udDataInds, 0]
        udResponsesNorth = udResponseVals[:, udDataInds, 1]
        udResponsesUp = udResponseVals[:, udDataInds, 2]
        responseUd.close()

        # Create design matrix.
        print("  Creating design matrix:")
        sys.stdout.flush()
        nE = self.numDataPoints
        nN = self.numDataPoints
        nU = self.numDataPoints
        self.design = numpy.zeros((self.numDesignRows, self.numImpulses),
                                  dtype=numpy.float64)
        self.design[0:nE,:] = numpy.transpose(self.llComp * llResponsesEast +
                                               self.udComp * udResponsesEast)
        self.design[nE:nE + nN,:] = numpy.transpose(
            self.llComp * llResponsesNorth + self.udComp * udResponsesNorth)
        self.design[nE + nN:nE + nN + nU,:] = numpy.transpose(
            self.llComp * llResponsesUp + self.udComp * udResponsesUp)

        return

    def matchCoords(self, coordsRef, coords):
        """Function to provide indices that match the given set of coordinates to a
        reference set.
        """

        """
    This is a lot easier if you have scipy.
    import scipy
    tree = scipy.spatial.cKDTree(coordsRef)
    (distances, inds) = tree.query(coords)
    """

        diff = coordsRef[:,:, None] - coords[:,:, None].transpose()
        dist = numpy.linalg.norm(diff, axis=1)
        inds = numpy.argmin(dist, axis=0)
        distances = dist[inds].diagonal()

        return (distances, inds)

    def readData(self):
        """Function to read data, coordinates, and standard deviations.
        """
        print("Reading data values:")
        sys.stdout.flush()

        coords = []
        data = []
        cov = []
        dataE = []
        dataN = []
        dataU = []
        covE = []
        covN = []
        covU = []
        self.dataNames = []

        f = open(self.dataFile, 'r')
        lines = f.readlines()
        self.numDataPoints = len(lines) - 1
        self.numDesignRows = 3 * self.numDataPoints
        for line in range(1, self.numDataPoints + 1):
            lineSplit = lines[line].split()
            x = float(lineSplit[1])
            y = float(lineSplit[2])
            z = float(lineSplit[3])
            self.dataNames.append(lineSplit[0])
            coords.append([x, y, z])
            vE = self.dataScale * float(lineSplit[4])
            vN = self.dataScale * float(lineSplit[5])
            vU = self.dataScale * float(lineSplit[6])
            dataE.append(vE)
            dataN.append(vN)
            dataU.append(vU)
            sigE = self.dataScale * float(lineSplit[7])
            sigN = self.dataScale * float(lineSplit[8])
            sigU = self.dataScale * float(lineSplit[9])
            covE.append(sigE * sigE)
            covN.append(sigN * sigN)
            covU.append(sigU * sigU)

        f.close()

        print("  Number of data locations: %i" % self.numDataPoints)
        print("  Number of rows in design matrix: %i" % self.numDesignRows)
        sys.stdout.flush()

        data = dataE + dataN + dataU
        cov = covE + covN + covU
        self.dataVals = numpy.array(data, dtype=numpy.float64)
        self.dataCov = numpy.array(cov, dtype=numpy.float64)
        self.dataCoords = numpy.array(coords, dtype=numpy.float64)

        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = SlipInvert()
    app.run()

# End of file
