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

## @file euler/euler

## @brief Python application to generate a spatial database based on
## specified Euler poles at a given set of points.

import os, re, sys

from pyre.applications.Script import Script as Application

class Euler(Application):
  """
  Python application to create dislocation/displacement BC for a
  specified Euler pole.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing Euler facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Euler facilities and properties.
    ##
    ## \b Properties
    ## @li \b spatial_dim Spatial dimension of coordinates.
    ## @li \b data_dim Dimension of data.
    ## @li \b points_file Filename of file containing point coordinates.
    ## *** Leave these out for now and replace with inventory of projector.
    ## @li \b points_datum Datum used for point coordinates.
    ## @li \b points_zone UTM zone used for point coordinates.
    ## *** Leave these out for now and replace with inventory of projector.
    ## @li \b points_offset_east East offset of point coords from UTM coords.
    ## @li \b points_offset_north North offset of point coords from UTM coords.
    ## @li \b euler_lat Latitude of Euler pole.
    ## @li \b euler_lon Longitude of Euler pole.
    ## @li \b euler_rot Rotation value for Euler pole (CCW positive).
    ##
    ## \b Facilities
    ## @li \b src_coordsys Coordinate system to convert from.
    ## @li \b dest_coord_sys Coordinate system to convert to.
    ## @li \b projector Projector for coordinate systems.
    ## @li None

    import pyre.inventory
    from pyre.units.angle import deg
    from pyre.units.length import m

    from spatialdata.geocoords.CSCart import CSCart
    srcCoordSys = pyre.inventory.facility("src_coordsys",
                                          family="src_coordsys", factory=CSCart)
    srcCoordSys.meta['tip'] = "Source coordinate system."

    from spatialdata.geocoords.CSGeo import CSGeo
    destCoordSys = pyre.inventory.facility("dest_coordsys",
                                          family="dest_coordsys", factory=CSGeo)
    destCoordSys.meta['tip'] = "Destination coordinate system."

    from spatialdata.geocoords.Projector import Projector
    projector = pyre.inventory.facility("projector",
                                        family="projector", factory=Projector)
    projector.meta['tip'] = "Coordinate system projector."

    # NOTE:  I don't think it makes sense to have a spatial dimension
    # other than 3 for this code, but I'm leaving it for now.
    spatialDim = pyre.inventory.int("spatial_dim", default=3)
    spatialDim.meta['tip'] = "Spatial dimension of coordinates."

    dataDim = pyre.inventory.int("data_dim", default=2)
    dataDim.meta['tip'] = "Dimension of data."

    bcType = pyre.inventory.str("bc_type", default="dislocation")
    bcType.meta['tip'] = "Type of BC (dislocation or displacement)."

    pointsFile = pyre.inventory.str("points_file", default="points.def")
    pointsFile.meta['tip'] = "Filename of file containing point coordinates."

    pointsSpatialDB = pyre.inventory.str("points_spatialdb",
                                         default="points.spatialdb")
    pointsSpatialDB.meta['tip'] = "Filename of output spatial database."

    # pointsDatum = pyre.inventory.str("points_datum", default="WGS84")
    # pointsDatum.meta['tip'] = "Datum used for point coordinates."

    # pointsZone = pyre.inventory.int("points_zone", default=11)
    # pointsZone.meta['tip'] = "UTM zone used for point coordinates."

    pointsOffsetEast = pyre.inventory.dimensional("points_offset_east",
                                                  default=0.0*m)
    pointsOffsetEast.meta['tip'] = \
                                 "East offset of point coords from UTM coords."

    pointsOffsetNorth = pyre.inventory.dimensional("points_offset_north",
                                                   default=0.0*m)
    pointsOffsetNorth.meta['tip'] = \
                                  "North offset of point coords from UTM coords."

    eulerLat = pyre.inventory.dimensional("euler_lat", default=0.0*deg)
    eulerLat.meta['tip'] = "Latitude of Euler pole."

    eulerLon = pyre.inventory.dimensional("euler_lon", default=0.0*deg)
    eulerLon.meta['tip'] = "Longitude of Euler pole."

    eulerRot = pyre.inventory.dimensional("euler_rot", default=0.0*deg)
    eulerRot.meta['tip'] = "Rotation of Euler pole (CCW positive)."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="euler"):
    Application.__init__(self, name)
    self.numPoints = 0
    self.pointsUTM = []
    self.normals = []
    return


  def main(self):
    self._readPoints()
    f = open(self.pointsSpatialDB, 'w')
    self._writeSpatialDBHead(f)
    self._genSpatialDB(f)
    f.close()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.srcCoordSys = self.inventory.srcCoordSys
    self.destCoordSys = self.inventory.destCoordSys
    self.projector = self.inventory.projector
    self.spatialDim = self.inventory.spatialDim
    self.dataDim = self.inventory.dataDim
    self.bcType = self.inventory.bcType
    self.pointsFile = self.inventory.pointsFile
    self.pointsSpatialDB = self.inventory.pointsSpatialDB
    # self.pointsDatum = self.inventory.pointsDatum
    # self.pointsZone = self.inventory.pointsZone
    self.pointsOffsetEast = self.inventory.pointsOffsetEast
    self.pointsOffsetNorth = self.inventory.pointsOffsetNorth
    self.eulerLat = self.inventory.eulerLat
    self.eulerLon = self.inventory.eulerLon
    self.eulerRot = self.inventory.eulerRot
    return


  def _writeSpatialDBHead(self, f):
    """
    Writes header portion of spatial database.
    """
    f.write('#SPATIAL.ascii 1\n')
    f.write('SimpleDB {\n')
    f.write('  num-values = 3\n')
    if self.bcType == 'dislocation':
      f.write('  value-names = left-lateral-slip   reverse-slip  fault-opening\n')
    else:
      f.write('  value-names = dof-0   dof-1  dof-2\n')
    f.write('  value-units = m m m\n')
    f.write('  num-locs = '+str(self.numPoints)+'\n')
    f.write('  data-dim = '+str(self.dataDim)+'\n')
    f.write('  space-dim = '+str(self.spatialDim)+'\n')
    f.write('  cs-data = cartesian {\n')
    f.write('    to-meters = 1.0\n')
    f.write('    space-dim = '+str(self.spatialDim)+'\n')
    f.write('  }\n')
    f.write('}\n')
    return


  def _genSpatialDB(self, f):
    """
    Computes dislocations/displacements from Euler pole and writes to spatial DB.
    """
    import numpy
    #  Need to fix from here.
    #  Create a numpy array from the pointsUTM list, but I need to save
    #  the original list.
    #  Use convert to convert to geographic, then loop over the points.
    #  First, compute displacement (velocity) vectors in global coordinates,
    #  then convert to (E, N, U) coordinates, and then, for dislocations,
    #  convert to fault-local coordinates and write results to spatial
    #  database.
    
    pointSize = 6
    coordUTM = [0.0, 0.0]
    normal = [0.0, 0.0, 0.0]
    for point in range(self.numPoints):
      cstart = point * pointSize
      nstart = cstart + 3
      coordUTM[0] = self.data[cstart] + self.pointsOffsetEast
      coordUTM[1] = self.data[cstart+1] + self.pointsOffsetNorth
      normal = self.data[nstart:nstart + 2]

    
  def _readPoints(self):
    f = file(self.pointsFile)
    for line in f.readlines():
      if not line.startswith('#'):
        data = line.split()
        # self.data.append([float(number) for number in line.split()])
        self.pointsUTM.append(float(data[0]+self.pointsOffsetEast)
        self.pointsUTM.append(float(data[1]+self.pointsOffsetNorth)
        self.normals.append(float(data[3:5])
        self.numPoints += 1
    f.close() 
    return
  
  
# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = Euler()
  app.run()

# End of file
