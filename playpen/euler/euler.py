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

import math
import numpy

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
    ## @li \b east_dir Vector defining east direction.
    ## @li \b north_dir Vector defining north direction.
    ## @li \b up_dir Vector defining up direction.
    ## @li \b euler_lat Latitude of Euler pole.
    ## @li \b euler_lon Longitude of Euler pole.
    ## @li \b euler_rot Rotation value for Euler pole (CCW positive).
    ##
    ## \b Facilities
    ## @li \b src_coordsys Coordinate system to convert from.
    ## @li \b dest_coordsys Coordinate system to convert to.
    ## @li \b projector Projector for coordinate systems.
    ## @li None

    import pyre.inventory
    from pyre.units.angle import deg
    from pyre.units.length import m

    from spatialdata.geocoords.CSGeoLocalCart import CSGeoLocalCart
    srcCoordSys = pyre.inventory.facility("src_coordsys",
                                          family="src_coordsys",
                                          factory=CSGeoLocalCart)
    srcCoordSys.meta['tip'] = "Source coordinate system."

    from spatialdata.geocoords.CSGeo import CSGeo
    destCoordSys = pyre.inventory.facility("dest_coordsys",
                                           family="dest_coordsys",
                                           factory=CSGeo)
    destCoordSys.meta['tip'] = "Destination coordinate system."

    from spatialdata.geocoords.Projector import Projector
    projector = pyre.inventory.facility("projector",
                                        family="projector", factory=Projector)
    projector.meta['tip'] = "Coordinate system projector."

    # NOTE:  I don't think it makes sense to have a spatial dimension
    # other than 3 for this code, but I'm leaving it for now.
    spatialDim = pyre.inventory.int("spatial_dim", default=3)
    spatialDim.meta['tip'] = "Spatial dimension of coordinates."

    dataDim = pyre.inventory.int("data_dim", default=3)
    dataDim.meta['tip'] = "Dimension of data."

    bcType = pyre.inventory.str("bc_type", default="dislocation")
    bcType.meta['tip'] = "Type of BC (dislocation or displacement)."

    pointsFile = pyre.inventory.str("points_file", default="points.def")
    pointsFile.meta['tip'] = "Filename of file containing point coordinates."

    pointsSpatialDB = pyre.inventory.str("points_spatialdb",
                                         default="points.spatialdb")
    pointsSpatialDB.meta['tip'] = "Filename of output spatial database."

    eastDir = pyre.inventory.list("east_dir", default=[1, 0, 0])
    eastDir.meta['tip'] = "East direction."

    northDir = pyre.inventory.list("north_dir", default=[0, 1, 0])
    northDir.meta['tip'] = "North direction."

    upDir = pyre.inventory.list("up_dir", default=[0, 0, 1])
    upDir.meta['tip'] = "Up direction."

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
    self.eulerPole = numpy.array([0.0, 0.0, 0.0])
    # Note that we use a mean radius since we are doing rotations on a
    # spherical Earth.
    self.earthRad = 63727954.77598
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
    self.eastDir = self.inventory.eastDir
    self.northDir = self.inventory.northDir
    self.upDir = self.inventory.upDir
    self.eulerLat = self.inventory.eulerLat
    self.eulerLon = self.inventory.eulerLon
    self.eulerRot = self.inventory.eulerRot

    lat = self.eulerLat.value
    lon = self.eulerLon.value
    rot = self.eulerRot.value
    clat = math.cos(lat)
    slat = math.sin(lat)
    clon = math.cos(lon)
    slon = math.sin(lon)
    # Note that the Euler pole already includes the earth radius term.
    self.eulerPole[0] = self.earthRad * rot * clat * clon
    self.eulerPole[1] = self.earthRad * rot * clat * slon
    self.eulerPole[2] = self.earthRad * rot * slat

    self.eastVec = numpy.array(self.eastDir)
    self.northVec = numpy.array(self.northDir)
    self.upVec = numpy.array(self.upDir)
    self.mesh2ENU = numpy.vstack((self.eastVec, self.northVec, self.upVec))
    
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

    # Get lat/lon values corresponding to UTM points.
    self.srcCoordSys.initialize()
    self.destCoordSys.initialize()

    from spatialdata.geocoords.Converter import convert
    pointsLL = numpy.array(self.pointsUTM).reshape(self.numPoints,
                                                   self.spatialDim-1)
    convert(pointsLL, self.destCoordSys, self.srcCoordSys)

    normalsArr = numpy.array(self.normals).reshape(self.numPoints,
                                                   self.spatialDim)
    
    iCount = 0
    for point in range(self.numPoints):
      velocity = self._euler2Velocity(pointsLL[point])
      if self.bcType == 'dislocation':
        self._localTrans(velocity, normalsArr[point])
      for comp in range(self.spatialDim):
        f.write(' %15e' % self.pointsUTM[iCount + comp])
      for comp in range(self.spatialDim):
        f.write(' %15e' % velocity[point, comp])
      f.write('\n')
      iCount += 3
    return


  def _localTrans(self, velocity, normalsArr):
    """
    Performs an in-place transformation on velocity vector to local coords.
    """
    # This function will need to partially duplicate the functionaliry of the
    # CellGeometry _orient2D function.
    mag = math.sqrt(numpy.dot(normalsArr, normalsArr))
    normalsArr /= mag

    # Along-strike direction -> cross product of up and normal
    alongStrike = numpy.cross(self.upVec, normalsArr)

    # Up-dip direction -> cross product of normal and along-strike
    upDip = numpy.cross(normalsArr, alongStrike)

    rot = numpy.vstack((alongStrike, upDip, normalsArr))
    velocity = numpy.dot(rot, velocity)
    return
    

  def _euler2Velocity(self, pointsLL):
    """
    Computes velocities in local Cartesian system from rotation about an
    Euler pole.
    """
    from pyre.units.angle import deg
    lonDeg = pointsLL[0]*deg
    latDeg = pointsLL[1]*deg
    lonPoint = lonDeg.value
    latPoint = latDeg.value
    clat = math.cos(latPoint)
    slat = math.sin(latPoint)
    clon = math.cos(lonPoint)
    slon = math.sin(lonPoint)
    pX = clat * clon
    pY = clat * slon
    pZ = slat
    pointPole = numpy.array([pX, pY, pZ])
    velGlobal = numpy.cross(self.eulerPole, pointPole)
    rotMatrix = self._makeRot(latPoint, lonPoint)
    velNED = numpy.dot(rotMatrix, velGlobal)
    # NOTE:  I should rearrange the rotation matrix so this transformation
    # isn't necessary.
    velENU = numpy.array([velNED(1), velNED(0), -velNED(2)])
    return velENU
    
      
  def _makeRot(self, latPoint, lonPoint):
    """
    Compute rotation matrix for a given geographic coordinates.
    """
    slat = math.sin(latPoint)
    clat = math.cos(latPoint)
    slon = math.sin(lonPoint)
    clon = math.cos(lonPoint)
    vec1 = [ -slat * clon, -slat * slon,  clat]
    vec2 = [        -slon,         clon,   0.0]
    vec3 = [ -clat * clon, -clat * slon, -slat]
    rotMatrix = numpy.array([vec1, vec2, vec3])
    return rotMatrix


  def _readPoints(self):
    """
    Reads point coordinates and normals from a file.
    """
    f = file(self.pointsFile)
    for line in f.readlines():
      if not line.startswith('#'):
        data = line.split()
        # self.data.append([float(number) for number in line.split()])
        self.pointsUTM.append(float(data[0:1]))
        self.normals.append(float(data[3:5]))
        self.numPoints += 1
    f.close() 
    return
  
  
# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = Euler()
  app.run()

# End of file
