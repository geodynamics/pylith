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

## @file greensfns/getcoords

## @brief Python application to get the coordinates for a set of vertex
## ID's and output them to a file. This is a very simple code that expects
## coordinates in Abaqus format as well as ID's in Abaqus nodeset format.

import math
import numpy
import os
import re
import glob
from pyre.units.time import s
from pyre.units.length import m

from pyre.applications.Script import Script as Application

class GetCoords(Application):
  """
  Python application to get the coordinates for a set of vertex
  ID's and output them to a file. This is a very simple code that expects
  coordinates in Abaqus format as well as ID's in Abaqus nodeset format.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing GetCoords facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing GetCoords facilities and properties.
    ##
    ## \b Properties
    ## @li \b coord_file  Name of file containing vertex coordinates.
    ## @li \b id_file  Name of file containing vertex ID's.
    ## @li \b output_file Name of output file with requested coordinates.
    ## \b Facilities
    ## @li None

    import pyre.inventory

    coordFile = pyre.inventory.str("coord_file", default="mesh.coords")
    coordFile.meta['tip'] = "Name of file containing vertex coordinates."

    idFile = pyre.inventory.str("id_file", default="mesh.ids")
    idFile.meta['tip'] = "Name of file containing vertex ID's."

    outputFile = pyre.inventory.str("output_file", default="mesh.ids")
    outputFile.meta['tip'] = "Name of output file with requested coordinates."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="getcoords"):
    Application.__init__(self, name)

    self.numIds = 0
    self.idList = []

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readIds()
    self._getCoords()
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
    self.coordFile = self.inventory.coordFile
    self.idFile = self.inventory.idFile
    self.outputFile = self.inventory.outputFile

    return
      

  def _readIds(self):
    """
    Function to read ID numbers.
    """
    f = open(self.idFile, 'r')
    lines = f.readlines()
    for line in lines:
      for entry in line.split(', '):
        self.idList.append(int(entry.rstrip(',\n')))

    self.numIds = len(self.idList)
    f.close()

    return


  def _getCoords(self):
    """
    Function to read a list of coordinates and output the coordinates if
    they correspond to one of the requested ID's.
    """
    f = open(self.coordFile, 'r')
    o = open(self.outputFile, 'w')
    lines = f.readlines()
    for line in lines:
      vertexId = int(line.split(', ')[0].rstrip(',\n'))
      for requestedId in self.idList:
        if vertexId == requestedId:
          x = line.split(', ')[1].rstrip(',')
          y = line.split(', ')[2].rstrip(',')
          z = line.split(', ')[3].rstrip(',\n')
          outLine = x + '   ' + y + '   ' + z + '\n'
          o.write(outLine)
          break

    f.close()
    o.close()

    return
  

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = GetCoords()
  app.run()

# End of file
