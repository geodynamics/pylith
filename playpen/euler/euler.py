#!/usr/bin/env nemesis
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

# @file euler/euler

# @brief Python application to generate a spatial database based on
# specified Euler poles at a given set of points.

import math
import numpy

from pythia.pyre.applications.Script import Script as Application


class Euler(Application):
    """Python application to create dislocation/displacement BC for a
    specified Euler pole.
    """

    class Inventory(Application.Inventory):
        """Python object for managing Euler facilities and properties.
        """

        # @class Inventory
        # Python object for managing Euler facilities and properties.
        ##
        # \b Properties
        # @li \b data_dim Dimension of data.
        # @li \b bc_type Boundary condition type (dislocation or displacement).
        # @li \b points_file Filename of file containing point coordinates.
        # @li \b points_spatialdb Filename of output spatial database.
        # @li \b bc_scale Scaling factor for output BC.
        # @li \b up_dir Vector defining up direction.
        # @li \b normal_dir General preferred normal direction.
        # @li \b euler_lat Latitude of Euler pole.
        # @li \b euler_lon Longitude of Euler pole.
        # @li \b euler_rot Rotation value for Euler pole (CCW positive).
        # @li \b dip_slip Allow dip-slip to accomodate non-strike-slip movement.
        # @li \b dip_cutoff Cutoff dip below which dip-slip movement is allowed.
        # @li \b x_min Minimum x-value for which to apply rotation.
        # @li \b x_max Maximum x-value for which to apply rotation.
        # @li \b y_min Minimum y-value for which to apply rotation.
        # @li \b y_max Maximum y-value for which to apply rotation.
        # @li \b z_min Minimum z-value for which to apply rotation.
        # @li \b z_max Maximum z-value for which to apply rotation.
        # @li \b default_values Values to use for out-of-range coordinates.
        ##
        # \b Facilities
        # @li \b src_coordsys Coordinate system to convert from.
        # @li \b dest_coordsys Coordinate system to convert to.

        import pythia.pyre.inventory
        from pythia.pyre.units.angle import deg
        from pythia.pyre.units.length import m
        from pythia.pyre.units.length import km

        from spatialdata.geocoords.CSGeoProj import CSGeoProj
        srcCoordSys = pythia.pyre.inventory.facility("src_coordsys",
                                                     family="src_coordsys",
                                                     factory=CSGeoProj)
        srcCoordSys.meta['tip'] = "Source coordinate system."

        from spatialdata.geocoords.CSGeo import CSGeo
        destCoordSys = pythia.pyre.inventory.facility("dest_coordsys",
                                                      family="dest_coordsys",
                                                      factory=CSGeo)
        destCoordSys.meta['tip'] = "Destination coordinate system."

        dataDim = pythia.pyre.inventory.int("data_dim", default=2)
        dataDim.meta['tip'] = "Dimension of data."

        bcType = pythia.pyre.inventory.str("bc_type", default="dislocation")
        bcType.meta['tip'] = "Type of BC (dislocation or displacement)."

        pointsFile = pythia.pyre.inventory.str(
            "points_file", default="points.def")
        pointsFile.meta['tip'] = "Filename of file containing point coordinates."

        pointsSpatialDB = pythia.pyre.inventory.str("points_spatialdb",
                                                    default="points.spatialdb")
        pointsSpatialDB.meta['tip'] = "Filename of output spatial database."

        bcScale = pythia.pyre.inventory.float("bc_scale", default=1.0)
        bcScale.meta['tip'] = "Scaling factor for output BC."

        upDir = pythia.pyre.inventory.list("up_dir", default=[0.0, 0.0, 1.0])
        upDir.meta['tip'] = "Up direction."

        normalDir = pythia.pyre.inventory.list(
            "normal_dir", default=[1.0, 0.0, 0.0])
        normalDir.meta['tip'] = "General preferred normal direction."

        eulerLat = pythia.pyre.inventory.dimensional(
            "euler_lat", default=0.0*deg)
        eulerLat.meta['tip'] = "Latitude of Euler pole."

        eulerLon = pythia.pyre.inventory.dimensional(
            "euler_lon", default=0.0*deg)
        eulerLon.meta['tip'] = "Longitude of Euler pole."

        eulerRot = pythia.pyre.inventory.dimensional(
            "euler_rot", default=0.0*deg)
        eulerRot.meta['tip'] = "Rotation of Euler pole (CCW positive)."

        dipSlip = pythia.pyre.inventory.bool("dip_slip", default=True)
        dipSlip.meta['tip'] = "Allow dip-slip to take up non-strike-slip movement."

        dipCutoff = pythia.pyre.inventory.dimensional(
            "dip_cutoff", default=75.0*deg)
        dipCutoff.meta['tip'] = "Cutoff dip below which dip-slip is allowed."

        xMin = pythia.pyre.inventory.dimensional("x_min", default=-1.0e8*m)
        xMin.meta['tip'] = "Minimum x-value for which to apply rotation."

        xMax = pythia.pyre.inventory.dimensional("x_max", default=1.0e8*m)
        xMax.meta['tip'] = "Maximum x-value for which to apply rotation."

        yMin = pythia.pyre.inventory.dimensional("y_min", default=-1.0e8*m)
        yMin.meta['tip'] = "Minimum y-value for which to apply rotation."

        yMax = pythia.pyre.inventory.dimensional("y_max", default=1.0e8*m)
        yMax.meta['tip'] = "Maximum y-value for which to apply rotation."

        zMin = pythia.pyre.inventory.dimensional("z_min", default=-1.0e8*m)
        zMin.meta['tip'] = "Minimum z-value for which to apply rotation."

        zMax = pythia.pyre.inventory.dimensional("z_max", default=1.0e8*m)
        zMax.meta['tip'] = "Maximum z-value for which to apply rotation."

        defaultValues = pythia.pyre.inventory.list("default_values",
                                                   default=[0.0, 0.0, 0.0])
        defaultValues.meta['tip'] = "Values used for out-of-range points."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="euler"):
        Application.__init__(self, name)
        self.numPoints = 0
        self.pointsUTM = []
        self.normals = []
        self.eulerPole = numpy.array([0.0, 0.0, 0.0], dtype=numpy.float64)
        # Note that we use a mean radius since we are doing rotations on a
        # spherical Earth.
        self.earthRad = 6372795.477598
        return

    def main(self):
        # import pdb
        # pdb.set_trace()
        self._readPoints()
        f = open(self.pointsSpatialDB, 'w')
        self._writeSpatialDBHead(f)
        self._genSpatialDB(f)
        f.close()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        Application._configure(self)
        self.srcCoordSys = self.inventory.srcCoordSys
        self.destCoordSys = self.inventory.destCoordSys
        self.dataDim = self.inventory.dataDim
        self.bcType = self.inventory.bcType
        self.pointsFile = self.inventory.pointsFile
        self.pointsSpatialDB = self.inventory.pointsSpatialDB
        self.bcScale = self.inventory.bcScale
        self.upDir = self.inventory.upDir
        self.normalDir = self.inventory.normalDir
        self.eulerLat = self.inventory.eulerLat
        self.eulerLon = self.inventory.eulerLon
        self.eulerRot = self.inventory.eulerRot
        self.spaceDim = self.srcCoordSys.spaceDim
        self.dipSlip = self.inventory.dipSlip
        self.dipCutoff = self.inventory.dipCutoff
        self.xMinVal = self.inventory.xMin.value
        self.xMaxVal = self.inventory.xMax.value
        self.yMinVal = self.inventory.yMin.value
        self.yMaxVal = self.inventory.yMax.value
        self.zMinVal = self.inventory.zMin.value
        self.zMaxVal = self.inventory.zMax.value
        self.defaultValues = self.inventory.defaultValues

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

        self.upVec = numpy.array([float(self.upDir[0]), float(self.upDir[1]),
                                  float(self.upDir[2])], dtype=numpy.float64)

        self.normalVec = numpy.array([float(self.normalDir[0]),
                                      float(self.normalDir[1]),
                                      float(self.normalDir[2])], dtype=numpy.float64)

        self.dipCutoffProj = abs(math.sin(self.dipCutoff.value))

        return

    def _writeSpatialDBHead(self, f):
        """Writes header portion of spatial database.
        """
        f.write('#SPATIAL.ascii 1\n')
        f.write('SimpleDB {\n')
        f.write('  num-values = 3\n')
        if self.bcType == 'dislocation':
            f.write(
                '  value-names = left-lateral-slip   reverse-slip  fault-opening\n')
        else:
            f.write('  value-names = dof-0   dof-1  dof-2\n')
        f.write('  value-units = m m m\n')
        f.write('  num-locs = '+str(self.numPoints)+'\n')
        f.write('  data-dim = '+str(self.dataDim)+'\n')
        f.write('  space-dim = '+str(self.spaceDim)+'\n')
        f.write('  cs-data = cartesian {\n')
        f.write('    to-meters = 1.0\n')
        f.write('    space-dim = '+str(self.spaceDim)+'\n')
        f.write('  }\n')
        f.write('}\n')
        return

    def _genSpatialDB(self, f):
        """Computes dislocations/displacements from Euler pole and writes to
        spatial DB.
        """

        # Get lat/lon values corresponding to UTM points.
        self.srcCoordSys.initialize()
        self.destCoordSys.initialize()

        from spatialdata.geocoords.Converter import convert
        pointsLL = numpy.array(self.pointsUTM,
                               dtype=numpy.float64).reshape(self.numPoints,
                                                            self.spaceDim)
        convert(pointsLL, self.destCoordSys, self.srcCoordSys)

        normalsArr = numpy.array(self.normals,
                                 dtype=numpy.float64).reshape(self.numPoints,
                                                              self.spaceDim)

        iCount = 0
        velocity = [0.0, 0.0, 0.0]
        for point in range(self.numPoints):
            inRange = self._testRange(self.pointsUTM[iCount:iCount+3])
            if inRange:
                velocity = self._euler2Velocity(pointsLL[point])
                if self.bcType == 'dislocation':
                    vlocal = self._localTrans(velocity, normalsArr[point])
                    velocity = self.bcScale * vlocal
            else:
                velocity[0] = self.bcScale * self.defaultValues[0]
                velocity[1] = self.bcScale * self.defaultValues[1]
                velocity[2] = self.bcScale * self.defaultValues[2]
            for dim in range(self.spaceDim):
                f.write(' %.12e' % self.pointsUTM[iCount + dim])
            for dim in range(self.spaceDim):
                f.write(' %.12e' % velocity[dim])
            f.write('\n')
            iCount += 3
        return

    def _testRange(self, point):
        """Checks to see if point is in range.
        """
        inRange = point[0] >= self.xMinVal and point[0] <= self.xMaxVal and \
            point[1] >= self.yMinVal and point[1] <= self.yMaxVal and \
            point[2] >= self.zMinVal and point[2] <= self.zMaxVal
        return inRange

    def _localTrans(self, velocity, normalsArr):
        """Performs a transformation on velocity vector to local coords.
        """
        # This function will need to partially duplicate the functionaliry of the
        # CellGeometry _orient2D function.
        mag = math.sqrt(numpy.dot(normalsArr, normalsArr))
        normalsArr /= mag

        normalTest = numpy.dot(normalsArr, self.normalVec)
        if normalTest < 0.0:
            normalsArr *= -1.0

        # Along-strike direction -> cross product of up and normal
        alongStrike = numpy.cross(self.upVec, normalsArr)
        mag = math.sqrt(numpy.dot(alongStrike, alongStrike))
        alongStrike /= mag

        # Up-dip direction -> cross product of normal and along-strike
        upDip = numpy.cross(normalsArr, alongStrike)
        mag = math.sqrt(numpy.dot(upDip, upDip))
        upDip /= mag

        rot = numpy.vstack((alongStrike, upDip, normalsArr))

        # Need to go through this section later to fix it for a generalized
        # coordinate setup.
        dip = numpy.dot(upDip, self.upVec)
        if self.dipSlip and abs(dip) <= self.dipCutoffProj:
            # Project slip onto strike-slip direction
            strikeSlipProj = numpy.dot(velocity, alongStrike)
            vstrikeSlip = strikeSlipProj * alongStrike

            # Horizontal normal movement is the difference between total velocity and
            # strike-slip velocity.
            vnormal = velocity - vstrikeSlip
            magHorizNormal = math.sqrt(
                vnormal[0]*vnormal[0]+vnormal[1]*vnormal[1])

            # Project horizontal normal movement onto dip-slip direction, then scale
            # so that horizontal components are equal to block-normal motion.
            dipSlipProj = numpy.dot(vnormal, upDip)
            vdipSlip = dipSlipProj * upDip
            magDipSlipHoriz = math.sqrt(vdipSlip[0]*vdipSlip[0] +
                                        vdipSlip[1]*vdipSlip[1])
            if magDipSlipHoriz > 0.0:
                multFac = magHorizNormal/magDipSlipHoriz
            else:
                multFac = 0.0
            vtotal = vstrikeSlip + multFac * vdipSlip
            vlocal = numpy.dot(rot, vtotal)
        else:
            vlocal = numpy.dot(rot, velocity)
        return vlocal

    def _euler2Velocity(self, pointsLL):
        """Computes velocities in local Cartesian system from rotation about an
        Euler pole.
        """
        from pythia.pyre.units.angle import deg
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
        pointPole = numpy.array([pX, pY, pZ], dtype=numpy.float64)
        velGlobal = numpy.cross(self.eulerPole, pointPole)
        rotMatrix = self._makeRot(latPoint, lonPoint)
        velNED = numpy.dot(rotMatrix, velGlobal)
        # NOTE:  I should rearrange the rotation matrix so this transformation
        # isn't necessary.
        velENU = numpy.array(
            [velNED[1], velNED[0], -velNED[2]], dtype=numpy.float64)
        return velENU

    def _makeRot(self, latPoint, lonPoint):
        """Compute rotation matrix for a given geographic coordinates.
        """
        slat = math.sin(latPoint)
        clat = math.cos(latPoint)
        slon = math.sin(lonPoint)
        clon = math.cos(lonPoint)
        vec1 = [-slat * clon, -slat * slon,  clat]
        vec2 = [-slon,         clon,   0.0]
        vec3 = [-clat * clon, -clat * slon, -slat]
        rotMatrix = numpy.array([vec1, vec2, vec3], dtype=numpy.float64)
        return rotMatrix

    def _readPoints(self):
        """Reads point coordinates and normals from a file.
        """
        f = file(self.pointsFile)
        for line in f.readlines():
            if not line.startswith('#'):
                data = line.split()
                # self.data.append([float(number) for number in line.split()])
                for dim in range(self.spaceDim):
                    self.pointsUTM.append(float(data[dim]))
                for dim in range(self.spaceDim):
                    self.normals.append(float(data[dim+self.spaceDim]))
                self.numPoints += 1
        f.close()
        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = Euler()
    app.run()

# End of file
