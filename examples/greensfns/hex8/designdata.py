#!/usr/bin/env python
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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file greensfns/designdata

## @brief Python application to create the data arrays necessary for an
## inversion, given data information and Green's function information.
## The matrices produced are the data matrix (vector), the data covariance
## matrix (diagonal at present), and the data design matrix (numObservations
## by numParameters).

import math
import numpy
import sys
from pyre.units.length import km
from pyre.units.length import m

from pyre.applications.Script import Script as Application

class DesignData(Application):
  """
  Python application to create the data arrays necessary for an
  inversion, given data information and Green's function information.
  The matrices produced are the data matrix (vector), the data covariance
  matrix (diagonal at present), and the data design matrix (numObservations
  by numParameters).
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing DesignData facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing DesignData facilities and properties.
    ##
    ## \b Properties
    ## @li \b data_input_file File containing data, locations, and stdDev.
    ## @li \b gf_metadata_file File containing metadata for GF.
    ## @li \b gfresponses_ll_root Root name for left-lateral GF responses.
    ## @li \b gfresponses_ud_root Root name for updip GF responses.
    ## @li \b data_output_file Output file for scaled data.
    ## @li \b cov_output_file Output file for scaled covariance matrix.
    ## @li \b design_output_file Output file for data design matrix.
    ## @li \b metadata_output_file Output file describing impulses and responses.
    ## @li \b data_scale Scaling factor to apply to data and stdDev.
    ## @li \b search_radius Radius from data center to search for GF.
    ## @li \b impulse_number_width Width of impulse number field.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    dataInputFile = pyre.inventory.str("data_input_file",
                                       default="data.txt")
    dataInputFile.meta['tip'] = "File containing data, locations, and stdDev."

    gfMetadataFile = pyre.inventory.str("gf_metadata_file",
                                        default="gf_metadata.txt")
    gfMetadataFile.meta['tip'] = "Name of file describing GF impulses."

    gfResponsesLlRoot = pyre.inventory.str("gfresponses_ll_root",
                                           default="gfresponse_ll.vtk")
    gfResponsesLlRoot.meta['tip'] = "Root name for left-lateral GF responses."

    gfResponsesUdRoot = pyre.inventory.str("gfresponses_ud_root",
                                           default="gfresponse_ud.vtk")
    gfResponsesUdRoot.meta['tip'] = "Root name for updip GF responses."

    dataOutputFile = pyre.inventory.str("data_output_file",
                                        default="data_vals.txt")
    dataOutputFile.meta['tip'] = "Output file for scaled data."

    covOutputFile = pyre.inventory.str("cov_output_file",
                                       default="data_cov.txt")
    covOutputFile.meta['tip'] = "Output file for scaled covriance matrix."

    designOutputFile = pyre.inventory.str("design_output_file",
                                          default="data_design.txt")
    designOutputFile.meta['tip'] = "Output file for data design matrix."

    metadataOutputFile = pyre.inventory.str("metadata_output_file",
                                            default="data_metadata.txt")
    metadataOutputFile.meta['tip'] = "Output file containing data metadata."

    dataScale = pyre.inventory.float("data_scale", default=1.0)
    dataScale.meta['tip'] = "Scaling factor to apply to data and stdDev."

    searchRadius = pyre.inventory.dimensional("search_radius", default=100.0*km)
    searchRadius.meta['tip'] = "Radius from data center to search for GF."

    impulseNumberWidth = pyre.inventory.int("impulse_number_width", default=5)
    impulseNumberWidth.meta['tip'] = "Width of impulse number field."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="designdata"):
    Application.__init__(self, name)

    self.numTotalImpulses = 0
    self.numUsedImpulses = 0
    self.numDataPoints = 0
    self.designRows = 0
    self.designColumns = 0
    self.numCells = 0
    self.numVertices = 0
    self.distanceScale = 0.0

    self.usedImpulses = []
    self.usedImpulsesLl = []
    self.usedImpulsesUd = []
    self.interpIndices = []
    self.dataNames = []
    self.dataCoords = None
    self.dataVals = None
    self.dataCov = None
    self.dataCenter = None
    self.design = None
    self.interpFuncs = None
    self.vertexCoords = None
    self.cellConnect = None

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readData()
    self._readMetadata()
    self._findImpulses()
    self._createInterp()
    self._writeMetadata()
    self._makeDesign()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    # import pdb
    # pdb.set_trace()

    # File info.
    self.dataInputFile = self.inventory.dataInputFile
    self.gfMetadataFile = self.inventory.gfMetadataFile
    self.gfResponsesLlRoot = self.inventory.gfResponsesLlRoot
    self.gfResponsesUdRoot = self.inventory.gfResponsesUdRoot
    self.dataOutputFile = self.inventory.dataOutputFile
    self.covOutputFile = self.inventory.covOutputFile
    self.designOutputFile = self.inventory.designOutputFile
    self.metadataOutputFile = self.inventory.metadataOutputFile

    # Data information
    self.dataScale = self.inventory.dataScale

    # Impulse information
    self.searchRadius = self.inventory.searchRadius.value
    self.impulseNumberWidth = self.inventory.impulseNumberWidth

    return
      

  def _readData(self):
    """
    Function to read data, coordinates, and standard deviations.
    """
    f = open(self.dataInputFile, 'r')
    lines = f.readlines()
    self.numDataPoints = len(lines) - 1
    self.designRows = 3 * self.numDataPoints
    coords = []
    data = []
    cov = []
    for line in range(1,self.numDataPoints + 1):
      lineSplit = lines[line].split()
      east = float(lineSplit[2])
      north = float(lineSplit[1])
      self.dataNames.append(lineSplit[0])
      coords.append([east, north])
      vE = self.dataScale * float(lineSplit[5])
      vN = self.dataScale * float(lineSplit[6])
      vU = self.dataScale * float(lineSplit[7])
      data.append(vE)
      data.append(vN)
      data.append(vU)
      sigE = self.dataScale * float(lineSplit[8])
      sigN = self.dataScale * float(lineSplit[9])
      sigU = self.dataScale * float(lineSplit[10])
      cov.append(sigE*sigE)
      cov.append(sigN*sigN)
      cov.append(sigU*sigU)

    f.close()

    print "Number of data points: %i" % self.numDataPoints
    print "Number of rows in design matrix: %i" % self.designRows
    sys.stdout.flush()
    self.dataVals = numpy.array(data, dtype=numpy.float64)
    self.dataCov = numpy.array(cov, dtype=numpy.float64)
    numpy.savetxt(self.dataOutputFile, self.dataVals)
    numpy.savetxt(self.covOutputFile, self.dataCov)
    self.dataCoords = numpy.array(coords, dtype=numpy.float64)
    self.dataCenter = numpy.mean(self.dataCoords, axis=0).reshape((1, 2))

    return


  def _readMetadata(self):
    """
    Function to read metadata.
    """
    self.metadata = numpy.loadtxt(self.gfMetadataFile, dtype=numpy.float64,
                                         skiprows=1)
    self.numTotalImpulses = self.metadata.shape[0]

    print "Total number of impulses:  %i" % self.numTotalImpulses
    sys.stdout.flush()

    return


  def _findImpulses(self):
    """
    Function to find impulses that lie within a given radius of the data center.
    """
    import scipy.spatial.distance
    
    impulseCoords = self.metadata[:,1:3]
    distance = scipy.spatial.distance.cdist(impulseCoords, self.dataCenter,
                                            'euclidean')

    for impulse in range(self.numTotalImpulses):
      if (distance[impulse] < self.searchRadius):
        self.usedImpulses.append(impulse)
        llFilename = self._getFilename(self.gfResponsesLlRoot, impulse)
        udFilename = self._getFilename(self.gfResponsesUdRoot, impulse)
        self.usedImpulsesLl.append(llFilename)
        self.usedImpulsesUd.append(udFilename)

    self.numUsedImpulses = len(self.usedImpulses)
    self.designColumns = 2 * self.numUsedImpulses
    print "Number of impulse locations used:  %i" % self.numUsedImpulses
    print "Number of columns in design matrix:  %i" % self.designColumns
    sys.stdout.flush()
    return

      
  def _getFilename(self, fileRoot, impulse):
    """
    Function to create a filename given the root filename and the impulse
    number.
    """
    impulseNum = int(impulse)
    impulseString = repr(impulseNum).rjust(self.impulseNumberWidth, '0')
    filename = fileRoot + "_t" + impulseString + ".vtk"

    return filename


  def _createInterp(self):
    """
    Function to find cell containing an observation point and create the
    corresponding interpolation functions.
    """
    from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
    from enthought.tvtk.api import tvtk

    # First read a VTK file to get geometry info
    reader = VTKFileReader()
    filename = self.usedImpulsesLl[0]
    reader.initialize(filename)
    data = reader.outputs[0]

    # Get cell and vertex info.
    cellVtk = data.get_cells()
    self.numCells = cellVtk.number_of_cells
    cellArray = cellVtk.to_array()
    self.vertexCoords = data._get_points().to_array()
    (self.numVertices, spaceDim) = self.vertexCoords.shape
    cellMatrix = cellArray.reshape(self.numCells,4)
    self.cellConnect = cellMatrix[:,1:4]

    meshRange = numpy.ptp(self.vertexCoords, axis=0)
    self.distanceScale = 0.1 * numpy.amax(meshRange)

    # Find cells enclosing each observation point.
    self.interpIndices = self._findCells()

    # Compute interpolation functions for each cell/observation point.
    self._computeInterp()

    return


  def _computeInterp(self):
    """
    Function to compute interpolation functions for a point lying in a
    triangular cell.
    """
    interpList = []
    for point in range(self.numDataPoints):
      cellNum = self.interpIndices[point]
      cellCoords = self.vertexCoords[self.cellConnect[cellNum,:],0:2]
      pointCoords = self.dataCoords[point,:]
      x1 = cellCoords[0,0]
      x2 = cellCoords[1,0]
      x3 = cellCoords[2,0]
      y1 = cellCoords[0,1]
      y2 = cellCoords[1,1]
      y3 = cellCoords[2,1]
      x = pointCoords[0]
      y = pointCoords[1]
      denom = x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1)
      l1 = x * (y3 - y2) + x2 * (y - y3) + x3 * (y2 - y)
      l2 = -(x * (y3 - y1) + x1 * (y - y3) + x3 * (y1 - y))
      l3 = x * (y2 - y1) + x1 * (y - y2) + x2 * (y1 - y)
      interpList.append([l1/denom, l2/denom, l3/denom])

    self.interpFuncs = numpy.array(interpList, dtype=numpy.float64)

    return


  def _findCells(self):
    """
    Function to find triangles enclosing a set of points.
    Brute force method for now.
    """

    indices = []
    for point in range(self.numDataPoints):
      pointCoords = self.dataCoords[point,:]
      cellFound = False
      for cell in range(self.numCells):
        cellCoords = self.vertexCoords[self.cellConnect[cell,:],0:2]
        inTriangle = self._inTriangle(pointCoords, cellCoords)
        if (inTriangle):
          indices.append(cell)
          cellFound = True
          break

      if (not cellFound):
        msg = 'Unable to find surrounding cell for data point # %i' % point
        raise ValueError(msg)

    return indices


  def _inTriangle(self, pointCoords, cellCoords):
    """
    Function to determine whether a point lies within a triangular cell.
    """

    v2 = pointCoords - cellCoords[0,:]

    if (math.sqrt(numpy.dot(v2, v2)) > self.distanceScale):
      return False

    v0 = cellCoords[2,:] - cellCoords[0,:]
    v1 = cellCoords[1,:] - cellCoords[0,:]

    dot00 = numpy.dot(v0, v0)
    dot01 = numpy.dot(v0, v1)
    dot02 = numpy.dot(v0, v2)
    dot11 = numpy.dot(v1, v1)
    dot12 = numpy.dot(v1, v2)

    invDenom = 1.0/(dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom

    return ((u > 0.0) and (v > 0.0) and (u + v < 1.0))
  
    
  def _makeDesign(self):
    """
    Function to create design matrix and write it to a file.
    """
    # Compute coefficients for each impulse/response pair.
    impulseList = []
    impulseTypes = [self.usedImpulsesLl, self.usedImpulsesUd]
    for impulseType in range(2):
      filelist = impulseTypes[impulseType]
      for impulse in range(self.numUsedImpulses):
        print 'Working on impulse # %i (%i)' % \
          (impulse, self.usedImpulses[impulse])
        sys.stdout.flush()
        impulseCoeffs = self._getCoeffs(filelist[impulse])
        impulseList.append(impulseCoeffs)
    self.design = numpy.transpose(numpy.array(impulseList, dtype=numpy.float64))
    numpy.savetxt(self.designOutputFile, self.design)

    return
    
      
  def _getCoeffs(self, vtkFile):
    """
    Function to compute all the coefficients from a particular impulse.
    """
    from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
    from enthought.tvtk.api import tvtk

    # Read VTK file
    reader = VTKFileReader()
    reader.initialize(vtkFile)
    data = reader.outputs[0]

    # Get computed displacements
    vertData = data._get_point_data()
    numVertDataArrays = vertData._get_number_of_arrays()
    displArray = None
    for arrayNum in range(numVertDataArrays):
      arrayName = vertData.get_array_name(arrayNum)
      if (arrayName == 'displacement'):
        displArray = vertData.get_array(arrayNum).to_array()
        break

    responseVals = []
    for dataPoint in range(self.numDataPoints):
      u, v, w = self._computeDispl(dataPoint, displArray)
      responseVals.append(u)
      responseVals.append(v)
      responseVals.append(w)

    return responseVals


  def _computeDispl(self, dataPoint, displArray):
    """
    Function to interpolate displacements to a given data point.
    """

    u = 0.0
    v = 0.0
    w = 0.0
    cellNum = self.interpIndices[dataPoint]
    for vertex in range(3):
      vertNum = self.cellConnect[cellNum, vertex]
      u += self.interpFuncs[dataPoint, vertex] * displArray[vertNum, 0]
      v += self.interpFuncs[dataPoint, vertex] * displArray[vertNum, 1]
      w += self.interpFuncs[dataPoint, vertex] * displArray[vertNum, 2]

    return (u, v, w)
  
    
  def _writeMetadata(self):
    """
    Function to write metadata describing impulses (columns) and data (rows).
    """
    f = open(self.metadataOutputFile, 'w')
    newLine = '\n'
    tab = '\t'
    
    # Write data information.
    dataDescr = 'Information for rows (data points), %i rows total\n' % self.designRows
    f.write(dataDescr)
    dataHead = 'Row #' + tab + 'Site-name' + tab + 'Component' + tab + \
      'X-coord' + tab + 'Y-coord' + newLine
    f.write(dataHead)
    components = ['X', 'Y', 'Z']
    dataFormat = '%i' + tab + '%s' + tab + '%s' + tab + '%e' + tab + '%e' + \
      newLine


    rowNum = 0
    for dataPoint in range(self.numDataPoints):
      coords = self.dataCoords[dataPoint, :]
      name = self.dataNames[dataPoint]
      for component in range(3):
        outLine = dataFormat % (rowNum, name, components[component],
                                coords[0], coords[1])
        f.write(outLine)
        rowNum += 1
    
    # Write impulse information.
    impulseDescr = newLine + newLine + \
                   'Information for columns (impulses), %i columns total\n' % self.designColumns
    f.write(impulseDescr)
    impulseHead = 'Column #' + tab + 'Impulse type' + tab + 'X-coord' + tab + \
                  'Y-coord' + tab + 'Z-coord' + tab + \
                  'Normal-X' + tab + 'Normal-Y' + tab + 'Normal-Z' + tab + \
                  'Strike-X' + tab + 'Strike-Y' + tab + 'Strike-Z' + tab + \
                  'Dip-X' + tab + 'Dip-Y' + tab + 'Dip-Z' + newLine
    f.write(impulseHead)
    impulseTypes = ['ll', 'ud']

    colNum = 0
    metadataOut = 12 * (tab + '%e')
    outFormat = '%i' + tab + '%s' + metadataOut + newLine
    for impulseType in range(2):
      impulseText = impulseTypes[impulseType]
      for impulse in self.usedImpulses:
        metadata = self.metadata[impulse, 1:13]
        outLine = outFormat % (colNum, impulseText,
                               metadata[0], metadata[1], metadata[2],
                               metadata[3], metadata[4], metadata[5],
                               metadata[6], metadata[7], metadata[8],
                               metadata[9], metadata[10], metadata[11])
        f.write(outLine)
        colNum += 1

    f.close()

    return
  

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = DesignData()
  app.run()

# End of file
