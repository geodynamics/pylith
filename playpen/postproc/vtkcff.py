#!/usr/bin/env python
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

# @file postproc/vtkcff

# @brief Python application to compute the Coulomb Failure Function difference
# for either a specified plane orientation or for optimally oriented planes.
# Differences are computed from a constant initial stress state, a specified
# state in a time series, or from the previous step in a time series.

import math
import numpy
import os
import re
import glob
from pythia.pyre.units.time import s

from pythia.pyre.applications.Script import Script as Application


class VtkCff(Application):
    """Python application to compute the Coulomb Failure Function (CFF) difference
    for either a specified plane orientation or for optimally oriented planes.
    Differences are computed from a zero initial stress state, a specified state
    in a time series, or from the previous step in a time series.
    """

    class Inventory(Application.Inventory):
        """Python object for managing VtkCff facilities and properties.
        """

        # @class Inventory
        # Python object for managing VtkCff facilities and properties.
        ##
        # \b Properties
        # @li \b stress_ref_mode Whether to reference stresses to a constant state, a selected state, or the previous state.
        # @li \b orientation_mode Compute CFF on predefined plane or optimally-oriented planes.
        # @li \b vtk_input_root Root filename for VTK input files.
        # @li \b vtk_output_root Root filename for VTK output files.
        # @li \b vtk_stress_index Index indicating which VTK field array contains stresses.
        # @li \b vtk_stress_components_order Indices corresponding to Sxx,Syy,Szz,Sxy,Syz,Sxz.
        # @li \b friction_coeff Coefficient of friction.
        # @li \b skempton_coeff Skempton's pore pressure coefficient B.
        # @li \b initial_state_index Initial state time step number for initial state mode.
        # @li \b constant_state_values Initial stress values to use with constant state mode.
        # @li \b cff_plane_normal Normal to plane on which to compute CFF
        # @li \b isotropic_poroelastic Use isotropic poroelastic model instead of constant apparent friction model.

        import pythia.pyre.inventory

        stressRefMode = pythia.pyre.inventory.str("stress_ref_mode",
                                           default="initial_state",
                                           validator=pythia.pyre.inventory.choice(["constant_state",
                                                                            "initial_state", "previous_state"]))
        stressRefMode.meta['tip'] = "Stress state against which to compute differences."

        orientationMode = pythia.pyre.inventory.str("orientation_mode",
                                             default="optimally_oriented",
                                             validator=pythia.pyre.inventory.choice(["optimally_oriented",
                                                                              "predefined_plane"]))
        orientationMode.meta['tip'] = "Compute CFF on predefined plane or optimally-oriented planes."

        vtkInputRoot = pythia.pyre.inventory.str("vtk_input_root",
                                          default="stress_t0001.vtk")
        vtkInputRoot.meta['tip'] = "Root filename for VTK input files."

        vtkOutputRoot = pythia.pyre.inventory.str("vtk_output_root", default="output.vtk")
        vtkOutputRoot.meta['tip'] = "Root filename for VTK output files."

        vtkStressIndex = pythia.pyre.inventory.int("vtk_stress_index", default=1)
        vtkStressIndex.meta['tip'] = "Index indicating which VTK field array contains stresses."

        vtkStressComponentsOrder = pythia.pyre.inventory.list("vtk_stress_components_order",
                                                       default=[0, 1, 2, 3, 4, 5])
        vtkStressComponentsOrder.meta['tip'] = "Indices corresponding to Sxx, Syy, Szz, Sxy, Syz, Sxz."

        frictionCoeff = pythia.pyre.inventory.float("friction_coeff", default=0.6)
        frictionCoeff.meta['tip'] = "Coefficient of friction."

        skemptonCoeff = pythia.pyre.inventory.float("skempton_coeff", default=0.5)
        skemptonCoeff.meta['tip'] = "Skempton's pore pressure coefficient B."

        initialStateIndex = pythia.pyre.inventory.int("initial_state_indes", default=0)
        initialStateIndex.meta['tip'] = "Initial state time step number for initial state mode."

        constantStateValues = pythia.pyre.inventory.list("constant_state_values",
                                                  default=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        constantStateValues.meta['tip'] = "Initial stress values to use with constant state mode."

        cffPlaneNormal = pythia.pyre.inventory.list("cff_plane_normal",
                                             default=[1.0, 0.0, 0.0])
        cffPlaneNormal.meta['tip'] = "Plane normal for predefined CFF plane."

        isotropicPoroelastic = pythia.pyre.inventory.bool("isotropic_poroelastic",
                                                   default=False)
        isotropicPoroelastic.meta['tip'] = "Use isotropic poroelastic model instead of constant apparent friction model."

    class TooFewFilesError(IOError):
        """Exception raised when not enough VTK input files are found.
        """

        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="vtkcff"):
        Application.__init__(self, name)
        self.vtkInputList = []
        self.numVtkInputFiles = 0
        self.vtkInputTimes = []
        self.timeStampWidth = 0
        self.refFileIndex = 0

        self.numVertsPerCell = 0
        self.numCells = 0
        self.cellsArray = numpy.array([0])
        self.verticesArray = numpy.array([0])
        self.numVerts = 0
        self.spaceDim = 0
        self.cellType = ""
        self.readMesh = False

        self.numStressPoints = 0
        return

    def main(self):
        import pdb
        pdb.set_trace()
        self._getFileInfo()
        self._cffLoop()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        Application._configure(self)
        import pdb
        pdb.set_trace()

        # Set up info for input files
        totalInputPath = os.path.normpath(
            os.path.join(os.getcwd(), self.inventory.vtkInputRoot))
        self.vtkInputDir = os.path.dirname(totalInputPath)
        baseInputName = os.path.basename(totalInputPath)
        baseInputNameLen = len(baseInputName)
        if baseInputName.endswith(".vtk"):
            baseInputNameStripped = baseInputName[0:baseInputNameLen - 4]
        else:
            baseInputNameStripped = baseInputName
        testFind = re.search('_t[0-9]*$', baseInputNameStripped)
        if testFind != None:
            timeInd = baseInputNameStripped.rfind(testFind.group(0))
            self.vtkInputRoot = baseInputNameStripped[0:timeInd]
        else:
            self.vtkInputRoot = baseInputNameStripped

        # Solution mode info
        self.stressRefMode = self.inventory.stressRefMode
        self.orientationMode = self.inventory.orientationMode
        self.initialStateIndex = self.inventory.initialStateIndex
        self.constantStateValues = self.inventory.constantStateValues
        self.isotropicPoroelastic = self.inventory.isotropicPoroelastic

        # Index information
        self.vtkStressIndex = self.inventory.vtkStressIndex
        self.vtkStressComponentsOrder = self.inventory.vtkStressComponentsOrder

        # Parameters
        self.frictionCoeff = self.inventory.frictionCoeff
        self.skemptonCoeff = self.inventory.skemptonCoeff
        self.cffPlaneNormal = numpy.array([float(self.inventory.cffPlaneNormal[0]),
                                           float(self.inventory.cffPlaneNormal[1]),
                                           float(self.inventory.cffPlaneNormal[2]),
                                           ], dtype=numpy.float64)

        # Set up info for output files
        totalOutputPath = os.path.normpath(os.path.join(
            os.getcwd(), self.inventory.vtkOutputRoot))
        self.vtkOutputDir = os.path.dirname(totalOutputPath)
        baseOutputName = os.path.basename(totalOutputPath)
        baseOutputNameLen = len(baseOutputName)
        if baseOutputName.endswith(".vtk"):
            baseOutputNameStripped = baseOutputName[0:baseOutputNameLen - 4]
        else:
            baseOutputNameStripped = baseOutputName
        testFind = re.search('_t[0-9]*$', baseOutputNameStripped)
        if testFind != None:
            timeInd = baseOutputNameStripped.rfind(testFind.group(0))
            self.vtkOutputRoot = baseOutputNameStripped[0:timeInd]
        else:
            self.vtkOutputRoot = baseOutputNameStripped

        return

    def _getFileInfo(self):
        """Find input files and set up filenames for input and output.
        """

        # Create list of input files and associated times
        fileString = self.vtkInputRoot + "_t[0-9]*.vtk"
        searchString = os.path.join(self.vtkInputDir, fileString)
        self.vtkInputList = glob.glob(searchString)
        self.numVtkInputFiles = len(self.vtkInputList)
        index1 = self.vtkInputList[0].rfind("_t")
        index2 = self.vtkInputList[0].rfind(".vtk")
        self.timeStampWidth = index2 - index1 - 2
        for vtkFile in self.vtkInputList:
            timeString = vtkFile[index1 + 2:index2]
            self.vtkInputTimes.append(float(timeString))

        # Determine time index from which to start processing
        if self.stressRefMode == "constant_state":
            self.startIndex = 0
        elif self.stressRefMode == "previous_state":
            self.startIndex = 1
        else:
            self.startIndex = self.initialStateIndex + 1

        # Raise exception if there aren't enough files to process
        numFilesToProcess = self.numVtkInputFiles - self.startIndex - 1
        if numFilesToProcess < 1:
            try:
                raise TooFewFilesError(numFilesToProcess)
            except TooFewFilesError as err:
                print('Not enough files found for search string:  ', searchString)
                print('Number of files found:  ', err.value)

        # Create output directory if it doesn't exist
        if not os.path.isdir(self.vtkOutputDir):
            os.mkdir(self.vtkOutputDir)

        return

    def _cffLoop(self):
        """Function to loop over input files, compute Cff relative to reference state,
        and write the results to output files.
        """

        # Define function to use, depending on whether we are computing CFF on a
        # predefined plane or an optimally-oriented plane.
        if self.orientationMode == "optimally_oriented":
            cffFunct = self._cffOop
        else:
            cffFunct = self._cffPredefined

        # Define reference stresses unless differences are computed between each
        # time step.
        if self.stressRefMode == "constant_state":
            testFile = os.path.join(self.vtkInputDir, self.vtkInputList[0])
            testStress = self._readStress(testFile)
            refStress = numpy.tile(self.constantStateValues,
                                   (self.numStressPoints, 1))
        elif self.stressRefMode == "initial_state":
            refFile = os.path.normpath(os.path.join(self.vtkInputDir,
                                                    self.vtkInputList[self.initialStateIndex]))
            refStress = self._readStress(refFile)

        # Loop over input VTK files.
        for fileInd in range(self.startIndex, self.numVtkInputFiles):
            timeStamp = self.vtkInputTimes[fileInd]
            timeStampOut = int(timeStamp)
            timeStampOutString = repr(timeStampOut).rjust(self.timeStampWidth, '0')
            outputFileName = self.vtkOutputRoot + "_t" + timeStampOutString + ".vtk"
            vtkOutputFile = os.path.join(self.vtkOutputDir, outputFileName)
            if self.stressRefMode == "previous_state":
                refFile = os.path.normpath(os.path.join(self.vtkInputDir,
                                                        self.vtkInputList[fileInd - 1]))
                refStress = self._readStress(refFile)
            newFile = os.path.normpath(os.path.join(self.vtkInputDir,
                                                    self.vtkInputList[fileInd]))
            newStress = self._readStress(newFile)
            cffFunct(refStress, newStress, vtkOutputFile)

        return

    def _readStress(self, vtkFile):
        """Function to read stresses from a file and store the info in a numpy array.
        """
        from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
        from enthought.tvtk.api import tvtk

        reader = VTKFileReader()
        reader.initialize(vtkFile)
        data = reader.outputs[0]

        # Get vertex and cell info if it hasn't already been done
        if not self.readMesh:
            cellVtk = data.get_cells()
            self.numVertsPerCell = cellVtk._get_max_cell_size()
            self.numCells = cellVtk.number_of_cells
            cellArray = cellVtk.to_array()
            self.cells = tvtk.CellArray()
            self.cells.set_cells(self.numCells, cellArray)
            self.vertArray = data._get_points().to_array()
            self.cellType = data.get_cell_type(0)
            (self.numVerts, self.spaceDim) = self.vertArray.shape
            self.readMesh = True

        # Get cell fields and extract stresses
        cellData = data._get_cell_data()
        numCellDataArrays = cellData._get_number_of_arrays()
        stress = cellData.get_array(self.stressIndex).to_array()
        (self.numStressPoints, numCols) = stress.shape

        sxx = stress[:, self.stressComponentsOrder[0]]
        syy = stress[:, self.stressComponentsOrder[1]]
        szz = stress[:, self.stressComponentsOrder[2]]
        sxy = stress[:, self.stressComponentsOrder[3]]
        syz = stress[:, self.stressComponentsOrder[4]]
        sxz = stress[:, self.stressComponentsOrder[5]]
        stressOrdered = numpy.column_stack((sxx, syy, szz, sxy, syz, sxz))

        return stressOrdered

    def _princStress(self, stress):
        """Function to compute 3D principal stresses and sort them.
        """
        stressMat = numpy.array([(stress[0], stress[3], stress[5]),
                                 (stress[3], stress[1], stress[4]),
                                 (stress[5], stress[4], stress[2])],
                                dtype=numpy.float64)
        (princStress, princAxes) = numpy.linalg.eigh(stressMat)
        idx = princStress.argsort()
        princStressOrdered = princStress[idx]
        princAxesOrdered = princAxes[:, idx]
        return princStressOrdered, princAxesOrdered

    def _CffOop(self, refStress, newStress, vtkOutputFile):
        """Function to compute CFF for optimally-oriented planes and output results to
        a file.
        """
        from enthought.tvtk.api import tvtk

        beta = math.atan(self.frictionCoeff) / 2.0
        sinBeta = math.sin(beta)
        cosBeta = math.cos(beta)
        sin2Beta = math.sin(2.0 * beta)
        cos2Beta = math.cos(2.0 * beta)
        sinCosBeta = sinBeta * cosBeta
        sinBetaSq = sinBeta * sinBeta
        cosBetaSq = cosBeta * cosBeta
        cff = numpy.empty((self.numStressPoints), dtype=numpy.float64)
        failDir1 = numpy.empty((self.numStressPoints, self.spaceDim),
                               dtype=numpy.float64)
        failDir2 = numpy.empty((self.numStressPoints, self.spaceDim),
                               dtype=numpy.float64)
        effFrictionCoeff = self.frictionCoeff * (1.0 - self.skemptonCoeff)
        presFac = -self.skemptonCoeff / 3.0
        vec1 = numpy.array([cosBeta, 0.0, sinBeta], dtype=numpy.float64)
        vec2 = numpy.array([cosBeta, 0.0, -sinBeta], dtype=numpy.float64)

        # Loop over stress points, compute total stress and the associated
        # principal stresses, as well as the stress difference, and then
        # compute CFF.
        for point in range(self.numStressPoints):
            totStress = newStress[point,:]
            deltaStress = newStress[point,:] - refStress[point,:]
            (totPrincStress, totPrincAxes) = self._princStress(totStress)
            # Rotate stress changes into principal axis coordinate system for
            # total stress.
            deltaStressMat = numpy.array(
                [(deltaStress[0], deltaStress[3], deltaStress[5]),
                 (deltaStress[3], deltaStress[1], deltaStress[4]),
                    (deltaStress[5], deltaStress[4], deltaStress[2])],
                dtype=numpy.float64)
            prod1 = numpy.dot(totPrincAxes, deltaStressMat)
            deltaStressRot = numpy.dot(prod1, totPrincAxes.transpose())
            s33 = deltaStressRot[0, 0] * sinBetaSq - \
                2.0 * deltaStressRot[0, 2] * sinCosBeta + \
                deltaStressRot[2, 2] * cosBetaSq
            s13 = 0.5 * (deltaStressRot[2, 2] - deltaStressRot[0, 0]) * sin2Beta + \
                deltaStressRot[0, 2] * cos2Beta
            deltaNormStress = deltaStress[0] + deltaStress[1] + deltaStress[2]
            # Compute CFF depending on which model is used
            if self.isotropicPoroelastic:
                cff[point] = s13 + self.frictionCoeff * \
                    (s33 + presFac * deltaNormStress)
            else:
                cff[point] = s13 + effFrictionCoeff * s33

            # Get failure planes by rotating vector in principal axis coordinate
            # system into global coordinate system.
            # NOTE:  make sure I'm applying the rotation the right way.
            failDir1[point,:] = numpy.dot(totPrincAxes.transpose(), vec1)
            failDir2[point,:] = numpy.dot(totPrincAxes.transpose(), vec2)

        # Set up mesh info for VTK file
        mesh = tvtk.UnstructuredGrid(points=self.vertArray)
        mesh.set_cells(self.cellType, self.cells)

        # Add output fields and write VTK file
        # For now, output cff as a scalar, failDir1 as a vector, and failDir2 as
        # a general array.
        cffName = "cff"
        failDir1Name = "failure_plane_1"
        failDir2Name = "failure_plane_2"
        mesh.cell_data.scalars = cff
        mesh.cell_data.scalars.name = cffName
        mesh.cell_data.vectors = failDir1
        mesh.cell_data.vectors.name = failDir1Name
        mesh.cell_data.add_array(failDir2)
        w = tvtk.UnstructuredGridWriter(file_name=vtkOutputFile, input=mesh)
        w.write()

        return

    def _CffPredefined(self, refStress, newStress, vtkOutputFile):
        """Function to compute CFF for predefined planes and output results to a file.
        """
        from enthought.tvtk.api import tvtk

        beta = math.atan(self.frictionCoeff) / 2.0
        sinBeta = math.sin(beta)
        cosBeta = math.cos(beta)
        sin2Beta = math.sin(2.0 * beta)
        cos2Beta = math.cos(2.0 * beta)
        sinCosBeta = sinBeta * cosBeta
        sinBetaSq = sinBeta * sinBeta
        cosBetaSq = cosBeta * cosBeta
        cff = numpy.empty((self.numStressPoints), dtype=numpy.float64)
        failDir1 = numpy.empty((self.numStressPoints, self.spaceDim),
                               dtype=numpy.float64)
        failDir2 = numpy.empty((self.numStressPoints, self.spaceDim),
                               dtype=numpy.float64)
        effFrictionCoeff = self.frictionCoeff * (1.0 - self.skemptonCoeff)
        presFac = -self.skemptonCoeff / 3.0
        vec1 = numpy.array([cosBeta, 0.0, sinBeta], dtype=numpy.float64)
        vec2 = numpy.array([cosBeta, 0.0, -sinBeta], dtype=numpy.float64)

        # Loop over stress points, compute total stress and the associated
        # principal stresses, as well as the stress difference, and then
        # compute CFF.
        for point in range(self.numStressPoints):
            totStress = newStress[point,:]
            deltaStress = newStress[point,:] - refStress[point,:]
            (totPrincStress, totPrincAxes) = self._princStress(totStress)
            # Rotate stress changes into principal axis coordinate system for
            # total stress.
            deltaStressMat = numpy.array(
                [(deltaStress[0], deltaStress[3], deltaStress[5]),
                 (deltaStress[3], deltaStress[1], deltaStress[4]),
                    (deltaStress[5], deltaStress[4], deltaStress[2])],
                dtype=numpy.float64)
            prod1 = numpy.dot(totPrincAxes, deltaStressMat)
            deltaStressRot = numpy.dot(prod1, totPrincAxes.transpose())
            s33 = deltaStressRot[0, 0] * sinBetaSq - \
                2.0 * deltaStressRot[0, 2] * sinCosBeta + \
                deltaStressRot[2, 2] * cosBetaSq
            s13 = 0.5 * (deltaStressRot[2, 2] - deltaStressRot[0, 0]) * sin2Beta + \
                deltaStressRot[0, 2] * cos2Beta
            deltaNormStress = deltaStress[0] + deltaStress[1] + deltaStress[2]
            # Compute CFF depending on which model is used
            if self.isotropicPoroelastic:
                cff[point] = s13 + self.frictionCoeff * \
                    (s33 + presFac * deltaNormStress)
            else:
                cff[point] = s13 + effFrictionCoeff * s33

            # Get failure planes by rotating vector in principal axis coordinate
            # system into global coordinate system.
            # NOTE:  make sure I'm applying the rotation the right way.
            failDir1[point,:] = numpy.dot(totPrincAxes.transpose(), vec1)
            failDir2[point,:] = numpy.dot(totPrincAxes.transpose(), vec2)

        # Set up mesh info for VTK file
        mesh = tvtk.UnstructuredGrid(points=self.vertArray)
        mesh.set_cells(self.cellType, self.cells)

        # Add output fields and write VTK file
        # For now, output cff as a scalar, failDir1 as a vector, and failDir2 as
        # a general array.
        cffName = "cff"
        failDir1Name = "failure_plane_1"
        failDir2Name = "failure_plane_2"
        mesh.cell_data.scalars = cff
        mesh.cell_data.scalars.name = cffName
        mesh.cell_data.vectors = failDir1
        mesh.cell_data.vectors.name = failDir1Name
        mesh.cell_data.add_array(failDir2)
        w = tvtk.UnstructuredGridWriter(file_name=vtkOutputFile, input=mesh)
        w.write()

        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = VtkCff()
    app.run()

# End of file
