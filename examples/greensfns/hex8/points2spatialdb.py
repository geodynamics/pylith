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

## @file greensfns/points2spatialdb

## @brief Python application to create a spatialdb using a given set of points.
## A subset of the points is read from a separate file, along with corresponding
## parameter values. A spatial database containing all points is created, with
## all values set to zero other than the subset.

import math
import numpy
from pyre.units.length import km
from pyre.units.length import m

from pyre.applications.Script import Script as Application

class Points2Spatialdb(Application):
  """
  Python application to create a spatialdb using a given set of points.
  A subset of the points is read from a separate file, along with corresponding
  parameter values. A spatial database containing all points is created, with
  all values set to zero other than the subset.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing Points2Spatialdb facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Points2Spatialdb facilities and properties.
    ##
    ## \b Properties
    ## @li \b point_coord_file File containing all coordinates.
    ## @li \b point_value_file File containing point subset with values.
    ## @li \b parameter_names List of parameter names.
    ## @li \b distance_tol Distance tolerance to determine coincident points.
    ##
    ## \b Facilities
    ## @li \b geometry  Geometry for output database.
    ## @li \b iohandler  Object for output database.

    import pyre.inventory

    pointCoordFile = pyre.inventory.str("point_coord_file",
                                        default="point_coords.txt")
    pointCoordFile.meta['tip'] = "File containing all coordinates."

    pointValueFile = pyre.inventory.str("point_value_file",
                                        default="point_values.txt")
    pointValueFile.meta['tip'] = "File containing point subset with values."

    parameterNames = pyre.inventory.list("parameter_names",
                                        default=['left-lateral-slip',
                                                 'reverse-slip',
                                                 'fault-opening'])
    parameterNames.meta['tip'] = "List of parameter names."

    parameterUnits = pyre.inventory.list("parameter_units",
                                        default=['m','m','m'])
    parameterUnits.meta['tip'] = "List of parameter units."

    distanceTol = pyre.inventory.dimensional("distance_tol",
                                             default=10.0*m)
    distanceTol.meta['tip'] = "Distance tolerance for coincident points."

    from spatialdata.spatialdb.generator.Geometry import Geometry
    geometry = pyre.inventory.facility("geometry", family="geometry",
                                       factory=Geometry)
    geometry.meta['tip'] = "Geometry for output database."

    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = pyre.inventory.facility("iohandler", family="simpledb_io",
                                        factory=SimpleIOAscii)
    iohandler.meta['tip'] = "Object for writing database."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="points2spatialdb"):
    Application.__init__(self, name)

    self.numTotalPoints = 0
    self.numSubsetPoints = 0
    self.numParameters = 0

    self.coincidentPoints = []

    self.totalPoints = None
    self.subsetPoints = None
    self.subsetValues = None

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readPoints()
    self._findCoincident()
    self._writeSpatialdb()
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
    self.pointCoordFile = self.inventory.pointCoordFile
    self.pointValueFile = self.inventory.pointValueFile

    # Parameter info.
    self.parameterNames = self.inventory.parameterNames
    self.parameterUnits = self.inventory.parameterUnits

    # Parameters.
    self.distanceTol = self.inventory.distanceTol.value

    # Spatialdb output facilities.
    self.geometry = self.inventory.geometry
    self.iohandler = self.inventory.iohandler

    return
      

  def _readPoints(self):
    """
    Function to read the file containing all points as well as the file
    containing the point subset with values.
    """
    self.totalPoints = numpy.loadtxt(self.pointCoordFile, dtype=numpy.float64)
    subsetData = numpy.loadtxt(self.pointValueFile, dtype=numpy.float64)
    self.subsetPoints = subsetData[:,0:3]
    self.subsetValues = subsetData[:,3:]
    self.numTotalPoints = self.totalPoints.shape[0]
    self.numSubsetPoints = self.subsetPoints.shape[0]
    self.numParameters = self.subsetValues.shape[1]
    if (self.numParameters != len(self.parameterNames)):
      msg = "Number of parameters not consistent with parameter names."
      raise ValueError(msg)

    return
        
    
  def _findCoincident(self):
    """
    Function to find points in the total set matching those in the subset.
    """
    import scipy.spatial.distance

    distance = scipy.spatial.distance.cdist(self.subsetPoints, self.totalPoints,
                                            'euclidean')
    minIndices = numpy.argmin(distance, axis=1)
    
    for point in range(self.numSubsetPoints):
      matchPoint = minIndices[point]
      if (distance[point, minIndices[point]] < self.distanceTol):
        self.coincidentPoints.append(matchPoint)
      else:
        msg = "No matching point found for subset point # %i." % point
        raise ValueError(msg)

    return


  def _writeSpatialdb(self):
    """
    Function to write out the spatial database, after inserting the desired
    values.
    """
    # Create info for each parameter.
    values = []
    for parameter in range(self.numParameters):
      paramVals = numpy.zeros(self.numTotalPoints, dtype=numpy.float64)
      paramVals[self.coincidentPoints] = self.subsetValues[:,parameter]
      paramInfo = {'name': self.parameterNames[parameter],
                   'units': self.parameterUnits[parameter],
                   'data': paramVals.flatten()}
      values.append(paramInfo)

    # Write database.
    data = {'points': self.totalPoints,
            'coordsys': self.geometry.coordsys,
            'data_dim': self.geometry.dataDim,
            'values': values}
    self.iohandler.write(data)

    return
  

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = Points2Spatialdb()
  app.run()

# End of file
