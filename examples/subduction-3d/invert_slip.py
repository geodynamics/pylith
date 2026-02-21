#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""Simple inversion script that uses PyLith-generated Green's functions.

This script is used as part of example step07. To perform the inversion you
must have already:
1.  Run example step06.
2.  Run the make_synthetic_gnssdisp.py script to generate synthetic data.
3.  Generated the step07 Green's functions (step07a and step07b).
Once you have performed the steps above, you can run this script. The
parameters are defined in invert_slip.cfg.

Usage:
./invert_slip.py
"""

import math
import numpy
import h5py
from pythia.pyre.units.angle import degree

from pythia.pyre.applications.Script import Script as Application


class SlipInvert(Application):
    """Python application to perform a linear inversion for slip using
    PyLith-generated Green's functions.
    """

    from pythia.pyre import inventory

    dataFile = inventory.str("data_file", default="data.txt")
    dataFile.meta['tip'] = "File with displ., locations, and stdDev."

    rake = inventory.dimensional("rake", default=90.0 * degree)
    rake.meta['tip'] = "Assumed rake angle."

    gfImpulsesLLFile = inventory.str("gfimpulses_ll_file", default="gfimpulse_ll.h5")
    gfImpulsesLLFile.meta['tip'] = "HDF5 file with left-lateral GF impulses."

    gfImpulsesRVFile = inventory.str("gfimpulses_rv_file", default="gfimpulse_ud.h5")
    gfImpulsesRVFile.meta['tip'] = "HDF5 file with updip GF impulses."

    gfResponsesLLFile = inventory.str("gfresponses_ll_file", default="gfresponse_ll.h5")
    gfResponsesLLFile.meta['tip'] = "HDF5 file with left-lateral GF responses."

    gfResponsesRVFile = inventory.str("gfresponses_ud_file", default="gfresponse_ud.h5")
    gfResponsesRVFile.meta['tip'] = "HDF5 file with updip GF responses."

    aPrioriValue = inventory.float("a_priori_value", default=0.0)
    aPrioriValue.meta['tip'] = "A priori value for parameters."

    penaltyWeightVals = inventory.list("penalty_weight_vals", default=[0.1, 0.5, 1.0, 5.0, 10.0])
    penaltyWeightVals.meta['tip'] = "List of penalty weights."

    dataScale = inventory.float("data_scale", default=1.0)
    dataScale.meta['tip'] = "Scaling factor to apply to data and stdDev."

    resultSummaryFile = inventory.str("result_summary_file", default='result_summary.txt')
    resultSummaryFile.meta['tip'] = "Text file summarizing inversion results."

    slipOutputFile = inventory.str("slip_output_file", default='predicted_slip.h5')
    slipOutputFile.meta['tip'] = "HDF5 file with predicted slip results."

    displacementOutputFile = inventory.str("displacement_output_file", default='predicted_displacement.h5')
    displacementOutputFile.meta['tip'] = "HDF5 file with predicted displacements."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="invert_slip"):
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

        self.summaryHead = "\n".join([
            "Penalty-weight",
            "Data-residual",
            "Weighted-data-residual",
            "Penalty-residual",
            "Weighted-penalty-residual",
            "Total-residual",
            "Total-weighted-residual",
        ])
        self.numSummaryCols = 7

        self.design = None

    def main(self):
        self.readData()
        self.readGreens()
        self.runInversions()

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        Application._configure(self)
        self.penaltyWeights = numpy.array(self.penaltyWeightVals, dtype=numpy.float64)
        self.numPenaltyWeights = self.penaltyWeights.shape[0]

        # Left-lateral and updip components from assumed rake.
        self.llComp = math.cos(self.rake.value)
        self.rvComp = math.sin(self.rake.value)

    def runInversions(self):
        """Function to run inversions using a range of penalty parameters.
        """
        print("Running inversions:")

        # Open output files.
        d = h5py.File(self.displacementOutputFile, 'w')
        s = h5py.File(self.slipOutputFile, 'w')

        # Write fault mesh and time info.
        summaryInfo = numpy.zeros((self.numPenaltyWeights, self.numSummaryCols), dtype=numpy.float64)
        cellDimF = 2
        timesF = self.penaltyWeights.reshape((self.numPenaltyWeights, 1, 1))
        timesF = s.create_dataset('time', data=timesF, maxshape=(None, 1, 1))
        topoF = s.create_dataset('viz/topology/cells', data=self.faultCells, dtype='d')
        topoF.attrs['cell_dim'] = numpy.int64(cellDimF)
        slipVec = numpy.zeros((self.numPenaltyWeights, self.numFaultVerts, 3), dtype=numpy.float64)
        slipAlongRake = numpy.zeros((self.numPenaltyWeights, self.numFaultVerts, 1), dtype=numpy.float64)

        # Write data mesh and time info.
        cellDimD = 0
        topolD = numpy.arange(self.numDataPoints, dtype=numpy.int64).reshape(self.numDataPoints, 1)
        timesD = self.penaltyWeights.reshape((self.numPenaltyWeights, 1, 1))
        timesD = d.create_dataset('time', data=timesD, maxshape=(None, 1, 1))
        topoD = d.create_dataset('viz/topology/cells', data=topolD, dtype='d')
        topoD.attrs['cell_dim'] = numpy.int64(cellDimD)

        predictedDisp = numpy.zeros((self.numPenaltyWeights, self.numDataPoints, 3), dtype=numpy.float64)

        # Rescale equations using data standard deviations.
        dataStdDev = numpy.sqrt(self.dataCov)
        dataStdDevInvDiag = numpy.diag(1.0 / dataStdDev)
        dataScaledDesign = numpy.dot(dataStdDevInvDiag, self.design)
        dataScaledVals = numpy.dot(dataStdDevInvDiag, self.dataVals)

        # Create a priori parameter vector.
        paramVec = self.aPrioriValue * numpy.ones(self.numImpulses, dtype=numpy.float64)

        # Regularization array is just the identity matrix.
        regArray = numpy.identity(self.numImpulses, dtype=numpy.float64)

        # Loop over inversions.
        for i_penalty, penWeight in self.penaltyWeights:
            print(f"  Working on inversion {i_penalty}, penalty weight = {penWeight:%g}")
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

            summaryInfo[i_penalty, 0] = penWeight
            summaryInfo[i_penalty, 1] = dataResidualNorm
            summaryInfo[i_penalty, 2] = dataWeightResidualNorm
            summaryInfo[i_penalty, 3] = penaltyResidualNorm
            summaryInfo[i_penalty, 4] = penaltyWeightResidualNorm
            summaryInfo[i_penalty, 5] = totalResidualNorm
            summaryInfo[i_penalty, 6] = totalWeightResidualNorm

            slipAlongRake[i_penalty, self.impulseInds, 0] = solution
            slipVec[i_penalty, self.impulseInds, 0] = self.llComp * solution
            slipVec[i_penalty, self.impulseInds, 1] = self.rvComp * solution
            predictedDisp[i_penalty,:,:] = predicted.reshape((self.numDataPoints, 3), order='F')

            print(f"    Data residual:              {dataResidualNorm:%e}")
            print(f"    Weighted data residual:     {dataWeightResidualNorm:%e}")
            print(f"    Penalty residual:           {penaltyResidualNorm:%e}")
            print(f"    Weighted penalty residual:  {penaltyWeightResidualNorm:%e}")
            print(f"    Total residual:             {totalResidualNorm:%e}")
            print(f"    Weighted total residual:    {totalWeightResidualNorm:%e}")

        numpy.savetxt(self.resultSummaryFile, summaryInfo, delimiter='\t', header=self.summaryHead)

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

    def readGreens(self):
        """Function to read impulse and response info from PyLith output files.
        """
        print("Reading Green's functions:")

        # Read impulses.
        print("  Reading left-lateral impulses:")
        impulsesLL = h5py.File(self.gfImpulsesLLFile, 'r')
        self.faultVertCoords = impulsesLL['geometry/vertices'][:]
        self.numFaultVerts = self.faultVertCoords.shape[0]
        self.faultCells = numpy.array(impulsesLL['viz/topology/cells'][:], dtype=numpy.int64)
        self.numFaultCells = self.faultCells.shape[0]
        llSlip = impulsesLL['vertex_fields/slip'][:,:, 1]
        llImpInds = numpy.nonzero(llSlip != 0.0)
        self.impulseCoords = self.faultVertCoords[llImpInds[1]]
        self.numImpulses = self.impulseCoords.shape[0]

        print(f"    Number of fault vertices:     {self.numFaultVerts}")
        print(f"    Number of impulses:           {self.numImpulses}")

        (distances, self.impulseInds) = self.matchCoords(self.faultVertCoords, self.impulseCoords)
        impulsesLL.close()

        print("  Reading updip impulses:")
        impulsesRv = h5py.File(self.gfImpulsesRVFile, 'r')
        rvCoords = impulsesRv['geometry/vertices'][:]
        rvSlip = impulsesRv['vertex_fields/slip'][:,:, 2]
        rvImpInds = numpy.nonzero(rvSlip != 0.0)
        rvCoordsUsed = rvCoords[rvImpInds[1]]
        (_, udCoordInds) = self.matchCoords(rvCoordsUsed, self.impulseCoords)
        rvCoordsUsed = rvCoordsUsed[udCoordInds,:]
        impulsesRv.close()

        # Read responses.
        print("  Reading left-lateral responses:")
        responseLL = h5py.File(self.gfResponsesLLFile, 'r')
        llResponseCoords = responseLL['geometry/vertices'][:]
        llResponseVals = responseLL['vertex_fields/displacement'][:]

        (_, llDataInds) = self.matchCoords(llResponseCoords, self.dataCoords)
        llResponsesEast = llResponseVals[:, llDataInds, 0]
        llResponsesNorth = llResponseVals[:, llDataInds, 1]
        llResponsesUp = llResponseVals[:, llDataInds, 2]
        responseLL.close()

        print("  Reading updip responses:")
        responseRv = h5py.File(self.gfResponsesRVFile, 'r')
        rvResponseCoords = responseRv['geometry/vertices'][:]
        responseRvVals = responseRv['vertex_fields/displacement'][:]
        rvResponseVals = responseRvVals[udCoordInds,:,:]

        (_, rvDataInds) = self.matchCoords(rvResponseCoords, self.dataCoords)
        rvResponsesEast = rvResponseVals[:, rvDataInds, 0]
        rvResponsesNorth = rvResponseVals[:, rvDataInds, 1]
        rvResponsesUp = rvResponseVals[:, rvDataInds, 2]
        responseRv.close()

        # Create design matrix.
        print("  Creating design matrix:")
        nE = self.numDataPoints
        nN = self.numDataPoints
        nU = self.numDataPoints
        self.design = numpy.zeros((self.numDesignRows, self.numImpulses), dtype=numpy.float64)
        self.design[0:nE,:] = numpy.transpose(self.llComp * llResponsesEast + self.rvComp * rvResponsesEast)
        self.design[nE:nE + nN,:] = numpy.transpose(self.llComp * llResponsesNorth + self.rvComp * rvResponsesNorth)
        self.design[nE + nN:nE + nN + nU,:] = numpy.transpose(self.llComp * llResponsesUp + self.rvComp * rvResponsesUp)

    def matchCoords(self, coordsRef, coords):
        """Function to provide indices that match the given set of coordinates to a
        reference set.

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

    def readDataOld(self):
        """Function to read data, coordinates, and standard deviations.
        """
        print("Reading data values:")

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

        f = open(self.dataFile, 'r', encoding="utf-8")
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

        print(f"  Number of data locations: {self.numDataPoints}")
        print(f"  Number of rows in design matrix: {self.numDesignRows}")

        data = dataE + dataN + dataU
        cov = covE + covN + covU
        self.dataVals = numpy.array(data, dtype=numpy.float64)
        self.dataCov = numpy.array(cov, dtype=numpy.float64)
        self.dataCoords = numpy.array(coords, dtype=numpy.float64)


    def readData(self):
        """Function to read data, coordinates, and standard deviations.
        """
        print("Reading data values:")

        data = numpy.loadtxt(self.dataFile, encoding="utf-8")
        self.dataNames = data[:,0]
        self.dataCoords = data[1:4,:]
        self.dataVals = self.dataScale * data[4:7]
        self.dataCov = self.dataScale * data[7:10]

        self.numDataPoints = data.shape[0]
        self.numDesignRows = 3 * self.numDataPoints

        print(f"  Number of data locations: {self.numDataPoints}")
        print(f"  Number of rows in design matrix: {self.numDesignRows}")


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = SlipInvert()
    app.run()
