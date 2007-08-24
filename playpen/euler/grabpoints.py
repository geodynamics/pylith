#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file grabpoints/grabpoints

## @brief Python application to grab a set of points specified in a pset
## file from a UCD file and write them to a file.

import math
import numpy

from pyre.applications.Script import Script as Application

class GrabPoints(Application):
  """
  Python application to grab a specified set of point coordinates and
  values from a UCD file.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing GrabPoints facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing GrabPoints facilities and properties.
    ##
    ## \b Properties
    ## @li \b pset_file Filename of file specifying vertex numbers.
    ## @li \b ucd_file Filename of input UCD file.
    ## @li \b point_output_file Filename of output set of points and normals.
    ## @li \b values_list List specifying position of desired attributes in UCD file.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    psetFile = pyre.inventory.str("pset_file", default="test.pset")
    psetFile.meta['tip'] = "Filename of pset file specifying vertex indices."

    ucdFile = pyre.inventory.str("ucd_file", default="test.inp")
    ucdFile.meta['tip'] = "Filename of ucd file containing mesh and attributes."

    pointOutputFile = pyre.inventory.str("point_output_file",
                                         default="points.coordnorm")
    pointOutputFile.meta['tip'] = "Filename of output coordinates and normals."

    valuesList = pyre.inventory.list("values_list", default=[1, 2, 3])
    valuesList.meta['tip'] = "Position of desired values in UCD attributes."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="grabpoints"):
    Application.__init__(self, name)
    self.numPoints = 0
    self.indices = []
    self.pointCoords = []
    return


  def main(self):
    import pdb
    pdb.set_trace()
    self._readPset()
    self._grabPoints()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.psetFile = self.inventory.psetFile
    self.ucdFile = self.inventory.ucdFile
    self.pointOutputFile = self.inventory.pointOutputFile
    self.valuesList = self.inventory.valuesList
    return


  def _readPset(self):
    """
    Reads vertex indices from a pset file.
    """
    f = file(self.psetFile)
    lines = f.readlines()
    line2 = lines[1]
    self.numPoints = int(line2.split()[2])
    numLines = len(lines)
    for lineCount in range(2, numLines):
      line = lines[lineCount]
      for number in line.split():
        self.indices.append(int(number))
    self.indices.sort()
    f.close() 
    return


  def _grabPoints(self):
    """
    Reads vertex indices from a pset file.
    """
    f = file(self.ucdFile)
    lines = f.readlines()
    fileLen = len(lines)
    firstline = lines[0].split()
    numVerts = int(firstline[0])
    numCells = int(firstline[1])
    numVertAttrs = int(firstline[2])
    vertInd = 0
    ucdInd = 1
    # Get vertex coordinates
    for lineCount in range(1, numVerts+1):
      vertex = self.indices[vertInd]
      if vertex == ucdInd:
        data = lines[lineCount].split()
	for dim in range(1,4):
          self.pointCoords.append(float(data[dim]))
        vertInd += 1
      ucdInd += 1

    # Skip elements and then start reading normals/values and write out
    # the selected values.
    o = open(self.pointOutputFile, 'w')
    lineBegin = 2 + numVerts + numCells + numVertAttrs
    lineEnd = lineBegin + numVerts
    vertInd = 0
    ucdInd = 1
    coordCount = 0
    normals = [0.0, 0.0, 0.0]
    v0 = int(self.valuesList[0])
    v1 = int(self.valuesList[1])
    v2 = int(self.valuesList[2])
    for lineCount in range(lineBegin, lineEnd):
      vertex = self.indices[vertInd]

      if vertex == ucdInd:
        data = lines[lineCount].split()
        normals = [float(data[v0]), float(data[v1]), float(data[v2])]

        for dim in range(3):
          o.write(' %15e' % self.pointCoords[coordCount + dim])

        for dim in range(3):
          o.write(' %15e' % normals[dim])

        o.write('\n')
        vertInd += 1
        coordCount += 3

      ucdInd += 1

    f.close() 
    o.close() 
    return
  
  
# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = GrabPoints()
  app.run()

# End of file
