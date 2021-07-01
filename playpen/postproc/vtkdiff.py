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

# @file postproc/vtkdiff

# @brief Python application to compute the difference between fields
# in a set of VTK files.

import math
import numpy
from pythia.pyre.units.time import s

from pythia.pyre.applications.Script import Script as Application


class VtkDiff(Application):
    """Python application to compute the difference between fields in a set of VTK
    files.
    """

    class Inventory(Application.Inventory):
        """Python object for managing VtkDiff facilities and properties.
        """

        # @class Inventory
        # Python object for managing VtkDiff facilities and properties.
        ##
        # \b Properties
        # @li \b vtk_input_root Root filename for VTK input files.
        # @li \b vtk_output_root Root filename for VTK output files.
        # @li \b scale_factor Scale factor to apply to results.

        import pythia.pyre.inventory

        vtkInputRoot = pythia.pyre.inventory.str("vtk_input_root", default="input.vtk")
        vtkInputRoot.meta['tip'] = "Root filename for VTK input files."

        vtkOutputRoot = pythia.pyre.inventory.str("vtk_output_root", default="output.vtk")
        vtkOutputRoot.meta['tip'] = "Root filename for VTK output files."

        scaleFactor = pythia.pyre.inventory.float("scale_factor", default=1.0)
        scaleFactor.meta['tip'] = "Scaling factor to apply to results."

    class TooFewFilesError(IOError):
        """Exception raised when not enough VTK input files are found.
        """

        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="vtkdiff"):
        Application.__init__(self, name)
        self.vtkInputList = []
        self.numVtkInputFiles = 0
        self.vtkInputTimes = []
        self.timeStampWidth = 0

        self.numVertsPerCell = 0
        self.numCells = 0
        self.numVerts = 0
        self.spaceDim = 0
        self.cellType = None
        self.readMesh = False
        return

    def main(self):
        # import pdb
        # pdb.set_trace()
        self._getFileInfo()
        self._computeDiffs()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        Application._configure(self)
        import os
        import re

        # Set up info for input files
        totalInputPath = os.path.normpath(os.path.join(os.getcwd(),
                                                       self.inventory.vtkInputRoot))
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

        self.scaleFactor = self.inventory.scaleFactor

        return

    def _getFileInfo(self):
        """Find input files and set up filenames for input and output.
        """
        import os
        import glob

        # Create list of input files and associated times
        fileString = self.vtkInputRoot + "_t[0-9]*.vtk"
        searchString = os.path.join(self.vtkInputDir, fileString)
        self.vtkInputList = glob.glob(searchString)
        self.vtkInputList.sort()
        self.numVtkInputFiles = len(self.vtkInputList)
        if self.numVtkInputFiles < 2:
            try:
                raise TooFewFilesError(self.numVtkInputFiles)
            except TooFewFilesError as err:
                print('Not enough files found for search string:  ', searchString)
                print('Number of files found:  ', err.value)
        index1 = self.vtkInputList[0].rfind("_t")
        index2 = self.vtkInputList[0].rfind(".vtk")
        self.timeStampWidth = index2 - index1 - 2
        for vtkFile in self.vtkInputList:
            timeString = vtkFile[index1 + 2:index2]
            self.vtkInputTimes.append(float(timeString))

        # Create output directory if it doesn't exist
        if not os.path.isdir(self.vtkOutputDir):
            os.mkdir(self.vtkOutputDir)

        return

    def _computeDiffs(self):
        """Function to loop over input files, compute differences in the fields, and
        write the results to output files.
        """
        import os

        # Loop over input VTK files.
        for fileInd in range(self.numVtkInputFiles - 1):
            timeStamp1 = self.vtkInputTimes[fileInd]
            timeStamp2 = self.vtkInputTimes[fileInd + 1]
            dt = (timeStamp2 - timeStamp1) / self.scaleFactor
            # We assume that differences correspond to midpoint of the 2 time steps.
            timeStampOut = int(0.5 * (timeStamp1 + timeStamp2))
            timeStampOutString = repr(timeStampOut).rjust(self.timeStampWidth, '0')
            #outputFileName = self.vtkOutputRoot + "_t" + timeStampOutString + ".vtk"
            outputFileName = self.vtkOutputRoot + "_t" + timeStampOutString + ".vtu"
            vtkOutputFile = os.path.join(self.vtkOutputDir, outputFileName)
            self._diffFiles(self.vtkInputList[fileInd],
                            self.vtkInputList[fileInd + 1],
                            vtkOutputFile,
                            dt)

        return

    def _diffFiles(self, vtkFile1, vtkFile2, vtkFileOut, dt):
        """Function to compute field differences between two VTK files, divide the
        differences by dt, and output the results to a new VTK file.
        """
        from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
        from enthought.tvtk.api import tvtk

        # Set up input files
        reader1 = VTKFileReader()
        reader2 = VTKFileReader()
        reader1.initialize(vtkFile1)
        reader2.initialize(vtkFile2)
        data1 = reader1.outputs[0]
        data2 = reader2.outputs[0]

        # Get vertex and cell info if it hasn't already been done
        if not self.readMesh:
            cellVtk = data1.get_cells()
            self.numVertsPerCell = cellVtk._get_max_cell_size()
            self.numCells = cellVtk.number_of_cells
            cellArray = cellVtk.to_array()
            self.cells = tvtk.CellArray()
            self.cells.set_cells(self.numCells, cellArray)
            self.vertArray = data1._get_points().to_array()
            self.cellType = data1.get_cell_type(0)
            (self.numVerts, self.spaceDim) = self.vertArray.shape
            self.readMesh = True

        # Set up mesh info for VTK file
        mesh = tvtk.UnstructuredGrid(points=self.vertArray)
        mesh.set_cells(self.cellType, self.cells)

        # Get vertex fields and compute differences if the fields exist
        vertData1 = data1._get_point_data()
        numVertDataArrays = vertData1._get_number_of_arrays()
        if numVertDataArrays != 0:
            vertData2 = data2._get_point_data()
            # This is very kludgy because I haven't yet figured out how to include
            # multiple scalar or vector fields, and I also don't know how to put in
            # a name for a general array (represented as a field).
            numScalarsUsed = 0
            numVectorsUsed = 0
            for vertDataArray in range(numVertDataArrays):
                array1 = vertData1.get_array(vertDataArray).to_array()
                (numPoints, numCols) = array1.shape
                arrayName = vertData1.get_array_name(vertDataArray) + "/dt"
                array2 = vertData2.get_array(vertDataArray).to_array()
                arrayOut = (array2 - array1) / dt
                # This is wrong if we have a scalar field with 3 components
                if (numCols == 3 and numVectorsUsed == 0):
                    mesh.point_data.vectors = arrayOut
                    mesh.point_data.vectors.name = arrayName
                    numVectorsUsed += 1
                elif numScalarsUsed == 0:
                    mesh.point_data.scalars = arrayOut
                    mesh.point_data.scalars.name = arrayName
                    numScalarsUsed += 1
                # Kludge to add a general array
                else:
                    mesh.point_data.add_array(arrayOut)

        # Get cell fields and compute differences if the fields exist
        cellData1 = data1._get_cell_data()
        numCellDataArrays = cellData1._get_number_of_arrays()
        if numCellDataArrays != 0:
            cellData2 = data2._get_cell_data()
            # This is very kludgy because I haven't yet figured out how to include
            # multiple scalar or vector fields, and I also don't know how to put in
            # a name for a general array (represented as a field).
            numScalarsUsed = 0
            numVectorsUsed = 0
            for cellDataArray in range(numCellDataArrays):
                array1 = cellData1.get_array(cellDataArray).to_array()
                (numPoints, numCols) = array1.shape
                arrayName = cellData1.get_array_name(cellDataArray) + "/dt"
                array2 = cellData2.get_array(cellDataArray).to_array()
                arrayOut = (array2 - array1) / dt
                # This is wrong if we have a scalar field with 3 components
                if (numCols == 3 and numVectorsUsed == 0):
                    mesh.cell_data.vectors = arrayOut
                    mesh.cell_data.vectors.name = arrayName
                    numVectorsUsed += 1
                elif numScalarsUsed == 0:
                    mesh.cell_data.scalars = arrayOut
                    mesh.cell_data.scalars.name = arrayName
                    numScalarsUsed += 1
                # Kludge to add a general array
                else:
                    mesh.cell_data.add_array(arrayOut)

        # Write results to VTK file
        #w = tvtk.UnstructuredGridWriter(file_name=vtkFileOut, input=mesh)
        w = tvtk.XMLDataSetWriter(file_name=vtkFileOut, input=mesh)
        w.write()

        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = VtkDiff()
    app.run()

# End of file
