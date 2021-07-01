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

# @file transform/transform

# @brief Python application to generate a spatial database based on
# a slip distribution for a region.  The slip distribution is converted
# back to global coordinates and then transformed back to local coordinates
# for a specified set of points.

import math
import numpy

from pythia.pyre.applications.Script import Script as Application


class Transform(Application):
    """Python application to create a slip distribution for a specified set of
    points.
    """

    class Inventory(Application.Inventory):
        """Python object for managing Transform facilities and properties.
        """

        # @class Inventory
        # Python object for managing Transform facilities and properties.
        ##
        # \b Properties
        # @li \b data_dim Dimension of data.
        # @li \b points_file Filename of file containing point coordinates.
        # @li \b segs_file Filename of file containing fault segments.
        # @li \b slip_scale Scaling factor for fault slip.
        # @li \b points_spatialdb Filename of output spatial database.
        # @li \b up_dir Vector defining up direction.
        # @li \b normal_dir General preferred normal direction.
        # @li \b dip_slip Allow dip-slip to accomodate non-strike-slip movement.
        # @li \b dip_cutoff Cutoff dip below which dip-slip movement is allowed.
        # @li \b x_min Minimum x-value for which to compute values.
        # @li \b x_max Maximum x-value for which to compute values.
        # @li \b y_min Minimum y-value for which to compute values.
        # @li \b y_max Maximum y-value for which to compute values.
        # @li \b z_min Minimum z-value for which to compute values.
        # @li \b z_max Maximum z-value for which to compute values.
        # @li \b default_values Values to use for out-of-range coordinates.
        ##
        # \b Facilities
        # @li \b src_coordsys Coordinate system to convert from.
        # @li \b dest_coordsys Coordinate system to convert to.

        import pythia.pyre.inventory
        from pythia.pyre.units.angle import deg
        from pythia.pyre.units.length import m
        from pythia.pyre.units.length import km

        dataDim = pythia.pyre.inventory.int("data_dim", default=2)
        dataDim.meta['tip'] = "Dimension of data."

        pointsFile = pythia.pyre.inventory.str(
            "points_file", default="points.def")
        pointsFile.meta['tip'] = "Filename of file containing point coordinates."

        segsFile = pythia.pyre.inventory.str("segs_file", default="segs.def")
        segsFile.meta['tip'] = "Filename of file containing fault segments."

        slipScale = pythia.pyre.inventory.float("slip_scale", default=1.0)
        slipScale.meta['tip'] = "Scaling factor for fault slip."

        pointsSpatialDB = pythia.pyre.inventory.str("points_spatialdb",
                                                    default="points.spatialdb")
        pointsSpatialDB.meta['tip'] = "Filename of output spatial database."

        upDir = pythia.pyre.inventory.list("up_dir", default=[0.0, 0.0, 1.0])
        upDir.meta['tip'] = "Up direction."

        normalDir = pythia.pyre.inventory.list(
            "normal_dir", default=[1.0, 0.0, 0.0])
        normalDir.meta['tip'] = "General preferred normal direction."

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

    def __init__(self, name="transform"):
        Application.__init__(self, name)
        self.numPoints = 0
        self.numSegs = 0
        self.segInfo = None
        self.pointsUTM = []
        self.normals = []
        return

    def main(self):
        # import pdb
        # pdb.set_trace()
        self._readPoints()
        self._readSegs()
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
        self.dataDim = self.inventory.dataDim
        self.pointsFile = self.inventory.pointsFile
        self.segsFile = self.inventory.segsFile
        self.slipScale = self.inventory.slipScale
        self.spaceDim = 3
        self.pointsSpatialDB = self.inventory.pointsSpatialDB
        self.upDir = self.inventory.upDir
        self.normalDir = self.inventory.normalDir
        self.dipSlip = self.inventory.dipSlip
        self.dipCutoff = self.inventory.dipCutoff
        self.xMinVal = self.inventory.xMin.value
        self.xMaxVal = self.inventory.xMax.value
        self.yMinVal = self.inventory.yMin.value
        self.yMaxVal = self.inventory.yMax.value
        self.zMinVal = self.inventory.zMin.value
        self.zMaxVal = self.inventory.zMax.value
        self.defaultValues = self.inventory.defaultValues

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
        f.write('  value-names = left-lateral-slip   reverse-slip  fault-opening\n')
        f.write('  value-units = m m m\n')
        f.write('  num-locs = '+str(self.numPoints)+'\n')
        f.write('  data-dim = '+str(self.dataDim)+'\n')
        f.write('  space-dim = '+str(self.spaceDim)+'\n')
        f.write('  cs-data = geo-projected {\n')
        f.write('    to-meters = 1.0\n')
        f.write('    ellipsoid = clrk66\n')
        f.write('    datum-horiz = NAD27\n')
        f.write('    datum-vert = mean sea level\n')
        f.write('    projector = projection {\n')
        f.write('      projection = utm\n')
        f.write('      units = m\n')
        f.write('      proj-options = +zone=11\n')
        f.write('    }\n')
        f.write('  }\n')
        f.write('}\n')
        return

    def _genSpatialDB(self, f):
        """Computes dislocations/displacements from Euler pole and writes to
        spatial DB.
        """

        normalsArr = numpy.array(self.normals,
                                 dtype=numpy.float64).reshape(self.numPoints,
                                                              self.spaceDim)

        points = numpy.array(self.pointsUTM, dtype=numpy.float64).reshape(self.numPoints,
                                                                          self.spaceDim)

        iCount = 0
        velocity = [0.0, 0.0, 0.0]
        for point in range(self.numPoints):
            inRange = self._testRange(points[point])
            if inRange:
                velocity = self._local2Global(points[point], normalsArr[point])
                vlocal = self._localTrans(velocity, normalsArr[point])
                velocity = vlocal
            else:
                velocity[0] = self.defaultValues[0]
                velocity[1] = self.defaultValues[1]
                velocity[2] = self.defaultValues[2]
            for dim in range(self.spaceDim):
                f.write(' %.12e' % points[point, dim])
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

    def _findSeg(self, point):
        """Finds segment to use with a given point.
        """
        px = point[0]
        py = point[1]
        """
    # This method does not seem to work too well.
    iSeg = -10
    for seg in range(self.numSegs):
      sx1 = self.segInfo[seg,0]
      sy1 = self.segInfo[seg,1]
      sx2 = self.segInfo[seg,2]
      sy2 = self.segInfo[seg,3]
      xmin = min(sx1,sx2)
      xmax = max(sx1,sx2)
      ymin = min(sy1,sy2)
      ymax = max(sy1,sy2)
      if sx2 == sx1:
        sSlope = 1.0e20
        pSlope = 0.0
      elif sy2 == sy1:
        sSlope = 0.0
        pSlope = 1.0e20
      else:
        sSlope = (sy2-sy1)/(sx2-sx1)
        pSlope = -1.0/sSlope
      pInt = py - pSlope * px
      sInt = sy1 - sSlope * sx1
      xInt = (pInt - sInt)/(sSlope - pSlope)
      yInt = (pInt*sSlope - sInt*pSlope)/(sSlope - pSlope)
      if xInt >= xmin and xInt <= xmax and yInt >= ymin and yInt >= ymax:
        iSeg = seg
        break
    """
        distMin = 1.0e20
        iSeg = -10
        for seg in range(self.numSegs):
            sx1 = self.segInfo[seg, 0]
            sy1 = self.segInfo[seg, 1]
            sx2 = self.segInfo[seg, 2]
            sy2 = self.segInfo[seg, 3]
            sxmid = 0.5*(sx1 + sx2)
            symid = 0.5*(sy1 + sy2)
            dx = px - sxmid
            dy = py - symid
            dist = math.sqrt(dx*dx+dy*dy)
            if dist < distMin:
                distMin = dist
                iSeg = seg
        return iSeg

    def _local2Global(self, points, normalsArr):
        """Converts initial local coordinate system to global.
        """
        iSeg = self._findSeg(points)
        sx1 = self.segInfo[iSeg, 0]
        sy1 = self.segInfo[iSeg, 1]
        sx2 = self.segInfo[iSeg, 2]
        sy2 = self.segInfo[iSeg, 3]
        dip = self.segInfo[iSeg, 4]
        ss = self.segInfo[iSeg, 5]
        ds = self.segInfo[iSeg, 6]
        ts = self.segInfo[iSeg, 7]
        dx = sx2 - sx1
        dy = sy2 - sy1
        asVec = numpy.array([dx, dy, 0.0], dtype=numpy.float64)
        mag = math.sqrt(numpy.dot(asVec, asVec))
        asVec /= mag
        horPerp = numpy.cross(self.upVec, asVec)
        mag = math.sqrt(numpy.dot(horPerp, horPerp))
        horPerp /= mag
        if dip > 90.0:
            horPerp *= -1.0
        testVec = numpy.cross(asVec, self.normalVec)
        if testVec[2] > 0.0:
            asVec *= -1.0
        dipCos = math.sin(math.radians(dip))
        if math.fabs(dipCos) != 1.0:
            r = math.sqrt(1.0/(1.0-dipCos*dipCos))
            udVec = numpy.array(
                [horPerp[0]/r, horPerp[1]/r, dipCos], dtype=numpy.float64)
        else:
            udVec = numpy.array([0.0, 0.0, dipCos], dtype=numpy.float64)
        normVec = numpy.cross(asVec, udVec)
        slipVec = numpy.array([ss, ds, ts], dtype=numpy.float64)
        rot1 = numpy.vstack((asVec, udVec, normVec))
        rot = rot1.transpose()
        velocity = numpy.dot(rot, slipVec)
        return velocity

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

    def _readSegs(self):
        """Reads fault segment info from a file.
        """
        segtmp = []
        f = file(self.segsFile)
        iCount = 1
        for line in f.readlines():
            if not line.startswith('#'):
                data = line.split()
                for dim in range(2):
                    segtmp.append(float(data[dim+2]))
                if iCount % 2 == 0:
                    segtmp.append(float(data[5]))
                    for dim in range(self.spaceDim):
                        segtmp.append(self.slipScale*float(data[dim+6]))
                    self.numSegs += 1
                iCount += 1
        f.close()
        self.segInfo = numpy.array(
            segtmp, dtype=numpy.float64).reshape(self.numSegs, 8)
        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = Transform()
    app.run()

# End of file
