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
    dataFile.meta["tip"] = "File with displ., locations, and stdDev."

    rake = inventory.dimensional("rake", default=90.0 * degree)
    rake.meta["tip"] = "Assumed rake angle."

    gfImpulsesLLFile = inventory.str("gfimpulses_ll_file", default="gfimpulse_ll.h5")
    gfImpulsesLLFile.meta["tip"] = "HDF5 file with left-lateral GF impulses."

    gfImpulsesRVFile = inventory.str("gfimpulses_rv_file", default="gfimpulse_ud.h5")
    gfImpulsesRVFile.meta["tip"] = "HDF5 file with updip GF impulses."

    gfResponsesLLFile = inventory.str("gfresponses_ll_file", default="gfresponse_ll.h5")
    gfResponsesLLFile.meta["tip"] = "HDF5 file with left-lateral GF responses."

    gfResponsesRVFile = inventory.str("gfresponses_rv_file", default="gfresponse_ud.h5")
    gfResponsesRVFile.meta["tip"] = "HDF5 file with updip GF responses."

    aPrioriValue = inventory.float("a_priori_value", default=0.0)
    aPrioriValue.meta["tip"] = "A priori value for parameters."

    penaltyWeightVals = inventory.list(
        "penalty_weight_vals", default=[0.1, 0.5, 1.0, 5.0, 10.0]
    )
    penaltyWeightVals.meta["tip"] = "List of penalty weights."

    dataScale = inventory.float("data_scale", default=1.0)
    dataScale.meta["tip"] = "Scaling factor to apply to data and stdDev."

    resultSummaryFile = inventory.str(
        "result_summary_file", default="result_summary.txt"
    )
    resultSummaryFile.meta["tip"] = "Text file summarizing inversion results."

    slipOutputFile = inventory.str("slip_output_file", default="predicted_slip.h5")
    slipOutputFile.meta["tip"] = "HDF5 file with predicted slip results."

    displacementOutputFile = inventory.str(
        "displacement_output_file", default="predicted_displacement.h5"
    )
    displacementOutputFile.meta["tip"] = "HDF5 file with predicted displacements."

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
        self.faultVertCoords = None
        self.faultCells = None

        self.numImpulses = 0
        self.impulseSlip = None

        self.summaryHead = "\t".join(
            [
                "Penalty-weight",
                "Data-residual",
                "Weighted-data-residual",
                "Penalty-residual",
                "Weighted-penalty-residual",
                "Total-residual",
                "Total-weighted-residual",
            ]
        )
        self.numSummaryCols = 7

        self.design = None

    def main(self):
        self.readData()
        self.readGreens()
        self.runInversions()

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory."""
        Application._configure(self)
        self.penaltyWeights = numpy.array(self.penaltyWeightVals, dtype=numpy.float64)
        self.numPenaltyWeights = self.penaltyWeights.shape[0]

        # Left-lateral and reverse components for assumed rake.
        self.llComp = math.cos(self.rake.value)
        self.rvComp = math.sin(self.rake.value)

    def runInversions(self):
        """Function to run inversions using a range of penalty parameters."""
        print("Running inversions:")

        # Open output files.
        h5Disp = h5py.File(self.displacementOutputFile, "w")
        h5Slip = h5py.File(self.slipOutputFile, "w")

        # Write fault mesh and time info.
        summaryInfo = numpy.zeros(
            (self.numPenaltyWeights, self.numSummaryCols), dtype=numpy.float64
        )
        cellDimF = 2
        timesF = self.penaltyWeights.reshape((self.numPenaltyWeights, 1, 1))
        timesF = h5Slip.create_dataset("time", data=timesF, maxshape=(None, 1, 1))
        vertsF = h5Slip.create_dataset("geometry/vertices", data=self.faultVertCoords)
        topoF = h5Slip.create_dataset(
            "viz/topology/cells", data=self.faultCells, dtype="d"
        )
        topoF.attrs["cell_dim"] = numpy.int64(cellDimF)
        slipVec = numpy.zeros(
            (self.numPenaltyWeights, self.numFaultVerts, 3), dtype=numpy.float64
        )
        slipAlongRake = numpy.zeros(
            (self.numPenaltyWeights, self.numFaultVerts, 1), dtype=numpy.float64
        )

        # Write data mesh and time info.
        cellDimD = 0
        topolD = numpy.arange(self.numDataPoints, dtype=numpy.int64).reshape(
            self.numDataPoints, 1
        )
        timesD = self.penaltyWeights.reshape((self.numPenaltyWeights, 1, 1))
        timesD = h5Disp.create_dataset("time", data=timesD, maxshape=(None, 1, 1))
        vertsD = h5Disp.create_dataset("geometry/vertices", data=self.dataCoords)
        topoD = h5Disp.create_dataset("viz/topology/cells", data=topolD, dtype="d")
        topoD.attrs["cell_dim"] = numpy.int64(cellDimD)

        predictedDisp = numpy.zeros(
            (self.numPenaltyWeights, self.numDataPoints, 3), dtype=numpy.float64
        )

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
        for i_penalty, penWeight in enumerate(self.penaltyWeights):
            print(f"  Working on inversion {i_penalty}, penalty weight = {penWeight:g}")
            paramScaledDesign = penWeight * regArray
            paramScaledData = penWeight * paramVec
            designMat = numpy.vstack((dataScaledDesign, paramScaledDesign))
            dataVec = numpy.hstack((dataScaledVals, paramScaledData))
            designMatTrans = numpy.transpose(designMat)
            genInv = numpy.dot(
                numpy.linalg.inv(numpy.dot(designMatTrans, designMat)), designMatTrans
            )
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

            slipAlongRake[i_penalty, :, 0] = numpy.dot(solution, self.impulseSlip)
            slipVec[i_penalty, :, 0] = self.llComp * numpy.dot(
                solution, self.impulseSlip
            )
            slipVec[i_penalty, :, 1] = self.rvComp * numpy.dot(
                solution, self.impulseSlip
            )
            predictedDisp[i_penalty, :, :] = predicted.reshape(
                (self.numDataPoints, 3), order="F"
            )

            print(f"    Data residual:              {dataResidualNorm:e}")
            print(f"    Weighted data residual:     {dataWeightResidualNorm:e}")
            print(f"    Penalty residual:           {penaltyResidualNorm:e}")
            print(f"    Weighted penalty residual:  {penaltyWeightResidualNorm:e}")
            print(f"    Total residual:             {totalResidualNorm:e}")
            print(f"    Weighted total residual:    {totalWeightResidualNorm:e}")

        numpy.savetxt(
            self.resultSummaryFile, summaryInfo, delimiter="\t", header=self.summaryHead
        )

        # Write results to HDF5 files.
        rakeSlip = h5Slip.create_dataset("vertex_fields/rake_slip", data=slipAlongRake)
        rakeSlip.attrs["vector_field_type"] = "scalar"
        vecSlip = h5Slip.create_dataset("vertex_fields/slip_vector", data=slipVec)
        vecSlip.attrs["vector_field_type"] = "vector"
        vecDisp = h5Disp.create_dataset("vertex_fields/disp_vec", data=predictedDisp)
        vecDisp.attrs["vector_field_type"] = "vector"

        h5Slip.close()
        h5Disp.close()

        from pylith.meshio.Xdmf import Xdmf

        xdmfWriter = Xdmf()
        xdmfWriter.write(self.slipOutputFile)
        xdmfWriter.write(self.displacementOutputFile)

    def readGreens(self):
        """Function to read impulse and response info from PyLith output files."""
        print("Reading Green's functions:")

        # Read impulses.
        print("  Reading left-lateral slip impulses:")
        h5LLImpulses = h5py.File(self.gfImpulsesLLFile, "r")
        self.faultVertCoords = h5LLImpulses["geometry/vertices"][:]
        self.numFaultVerts = self.faultVertCoords.shape[0]
        self.faultCells = numpy.array(
            h5LLImpulses["viz/topology/cells"][:], dtype=numpy.int64
        )
        llSlip = h5LLImpulses["vertex_fields/slip"][:, :, 1]
        self.numImpulses = llSlip.shape[0]
        h5LLImpulses.close()
        print(f"    Number of fault vertices:     {self.numFaultVerts}")
        print(f"    Number of impulses:           {self.numImpulses}")

        print("  Reading reverse slip impulses:")
        h5RvImpulses = h5py.File(self.gfImpulsesRVFile, "r")
        rvSlip = h5RvImpulses["vertex_fields/slip"][:, :, 2]
        h5RvImpulses.close()

        # Read responses.
        print("  Reading responses for left-lateral slip impulses:")
        h5LLResponse = h5py.File(self.gfResponsesLLFile, "r")
        llResponseValues = h5LLResponse["vertex_fields/displacement"][:]

        llResponsesEast = llResponseValues[:, :, 0]
        llResponsesNorth = llResponseValues[:, :, 1]
        llResponsesUp = llResponseValues[:, :, 2]
        h5LLResponse.close()

        print("  Reading responses for reverse slip impulses:")
        h5RvResponse = h5py.File(self.gfResponsesRVFile, "r")
        rvResponseValues = h5RvResponse["vertex_fields/displacement"][:]
        rvResponsesEast = rvResponseValues[:, :, 0]
        rvResponsesNorth = rvResponseValues[:, :, 1]
        rvResponsesUp = rvResponseValues[:, :, 2]
        h5RvResponse.close()

        # Slip impulse in rake direction
        self.impulseSlip = self.llComp * llSlip + self.rvComp * rvSlip

        # Create design matrix.
        print("  Creating design matrix:")
        nE = self.numDataPoints
        nN = self.numDataPoints
        nU = self.numDataPoints
        self.design = numpy.zeros(
            (self.numDesignRows, self.numImpulses), dtype=numpy.float64
        )
        self.design[0:nE, :] = numpy.transpose(
            self.llComp * llResponsesEast + self.rvComp * rvResponsesEast
        )
        self.design[nE : nE + nN, :] = numpy.transpose(
            self.llComp * llResponsesNorth + self.rvComp * rvResponsesNorth
        )
        self.design[nE + nN : nE + nN + nU, :] = numpy.transpose(
            self.llComp * llResponsesUp + self.rvComp * rvResponsesUp
        )

    def matchCoords(self, coordsRef, coords):
        """Function to provide indices that match the given set of coordinates to a
        reference set.

        This is a lot easier if you have scipy.
        import scipy
        tree = scipy.spatial.cKDTree(coordsRef)
        (distances, inds) = tree.query(coords)
        """
        diff = coordsRef[:, :, None] - coords[:, :, None].T
        dist = numpy.linalg.norm(diff, axis=1)
        indices = numpy.argmin(dist, axis=0)
        distances = dist[indices].diagonal()
        return (distances, indices)

    def readData(self):
        """Function to read data, coordinates, and standard deviations."""
        print("Reading data values:")

        data = numpy.genfromtxt(
            self.dataFile,
            encoding="utf-8",
            delimiter=",",
            names=True,
            dtype=None,
        )
        self.dataNames = data["Station"]
        self.dataCoords = numpy.stack((data["X"], data["Y"], data["Z"])).T
        self.dataVals = (
            self.dataScale * numpy.stack((data["UEast"], data["UNorth"], data["UUp"]))
        ).ravel()
        self.dataCov = (
            self.dataScale
            * numpy.stack(
                (data["SigEast"] ** 2, data["SigNorth"] ** 2, data["SigUp"] ** 2)
            )
        ).ravel()

        self.numDataPoints = data.shape[0]
        self.numDesignRows = 3 * self.numDataPoints

        print(f"  Number of data locations: {self.numDataPoints}")
        print(f"  Number of rows in design matrix: {self.numDesignRows}")


# ----------------------------------------------------------------------
if __name__ == "__main__":
    app = SlipInvert()
    app.run()
