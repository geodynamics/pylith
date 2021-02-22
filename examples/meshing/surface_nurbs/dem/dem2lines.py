#!/usr/bin/env nemesis

# @file dem2lines.py

# @brief Python application to read an ASCII (x,y,z) representation of a DEM
# and create a set of lines suitable for use in Cubit. The original DEM is
# resampled to a coarser resolution outside a specified region.

import math
import numpy
import pdb

from pythia.pyre.applications.Script import Script as Application


class Dem2Lines(Application):
    """
    Python application to read an ASCII (x,y,z) representation and
    create a set of lines suitable for use in Cubit. The original DEM is
    resampled to a coarser resolution outside a specified region.
    The DEM is assumed to be ordered by rows (left to right).
    """

    class Inventory(Application.Inventory):
        """
        Python object for managing Dem2Lines facilities and properties.
        """

        # @class Inventory
        # Python object for managing Dem2Lines facilities and properties.
        ##
        # \b Properties
        # @li \b input_dem Input DEM file (ASCII).
        # @li \b vtk_output_file Filename of VTK output file.
        # @li \b master_journal Output master journal file.
        # @li \b u_line_prefix Prefix for u (east) output lines.
        # @li \b u_line_journal Output journal file for u (east) output lines.
        # @li \b v_line_prefix Prefix for v (north) output lines.
        # @li \b v_line_journal Output journal file for v (north) output lines.
        # @li \b acis_filename Name of ACIS output file created by Cubit.
        # @li \b x_min Minimum x-value for full resolution.
        # @li \b x_max Maximum x-value for full resolution.
        # @li \b y_min Minimum y-value for full resolution.
        # @li \b y_max Maximum y-value for full resolution.
        # @li \b skip_interval Increment by which to increase points skipped.

        import pythia.pyre.inventory
        from pythia.pyre.units.angle import degree

        inputDem = pythia.pyre.inventory.str("input_dem", default="DEM.txt")
        inputDem.meta['tip'] = "Input DEM file."

        vtkOutputFile = pythia.pyre.inventory.str("vtk_output_file", default="DEM.vtk")
        vtkOutputFile.meta['tip'] = "VTK output file."

        masterJournal = pythia.pyre.inventory.str("master_journal", default="mktopo.jou")
        masterJournal.meta['tip'] = "Filename of output master journal file."

        uLinePrefix = pythia.pyre.inventory.str("u_line_prefix", default="u_lines")
        uLinePrefix.meta['tip'] = "Prefix for u (east) lines."

        uLineJournal = pythia.pyre.inventory.str("u_line_journal", default="u_lines.jou")
        uLineJournal.meta['tip'] = "Outputjournal file for u (east) lines."

        vLinePrefix = pythia.pyre.inventory.str("v_line_prefix", default="v_lines")
        vLinePrefix.meta['tip'] = "Prefix for v (north) lines."

        vLineJournal = pythia.pyre.inventory.str("v_line_journal", default="v_lines.jou")
        vLineJournal.meta['tip'] = "Outputjournal file for v (north) lines."

        acisFilename = pythia.pyre.inventory.str("acis_filename", default="topo.sab")
        acisFilename.meta['tip'] = "Name of ACIS output file created by Cubit."

        xMin = pythia.pyre.inventory.float("x_min", default=-1.0)
        xMin.meta['tip'] = "Minimum x-value for full resolution."

        xMax = pythia.pyre.inventory.float("x_max", default=1.0)
        xMax.meta['tip'] = "Maximum x-value for full resolution."

        yMin = pythia.pyre.inventory.float("y_min", default=-1.0)
        yMin.meta['tip'] = "Minimum y-value for full resolution."

        yMax = pythia.pyre.inventory.float("y_max", default=1.0)
        yMax.meta['tip'] = "Maximum y-value for full resolution."

        skipInterval = pythia.pyre.inventory.int("skip_interval", default=1)
        skipInterval.meta['tip'] = "Increment by which to increase points skipped."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="dem2lines"):
        Application.__init__(self, name)
        self.numXIn = 0
        self.numYIn = 0
        self.numZIn = 0
        self.numXOut = 0
        self.numYOut = 0
        self.numZOut = 0
        self.xIn = None
        self.yIn = None
        self.zIn = None
        self.xOut = None
        self.yOut = None
        self.zOut = None

        return

    def main(self):
        # pdb.set_trace()
        self._readDem()
        self._resampleDem()
        self._writeCubitJournals()
        self._writeDemVtk()

        return

    # PRIVATE METHODS /////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        Application._configure(self)

        # Filenames
        self.inputDem = self.inventory.inputDem
        self.vtkOutputFile = self.inventory.vtkOutputFile
        self.masterJournal = self.inventory.masterJournal
        self.uLinePrefix = self.inventory.uLinePrefix
        self.uLineJournal = self.inventory.uLineJournal
        self.vLinePrefix = self.inventory.vLinePrefix
        self.vLineJournal = self.inventory.vLineJournal
        self.acisFilename = self.inventory.acisFilename

        # Parameters
        self.xMin = self.inventory.xMin
        self.xMax = self.inventory.xMax
        self.yMin = self.inventory.yMin
        self.yMax = self.inventory.yMax
        self.skipInterval = self.inventory.skipInterval

        return

    def _readDem(self):
        """
        Read coordinates defining DEM and create vectors of x, y, and z values.
        """

        # Load each coordinate as a numpy array.
        x, y, z = numpy.loadtxt(self.inputDem, dtype=numpy.float64, unpack=True)

        self.numZIn = len(z)
        if (y[0] == y[1]):
            # Ordered by rows.
            self.numXIn = max(numpy.argmax(x) + 1, numpy.argmin(x) + 1)
            self.xIn = x[0:self.numXIn]
            self.numYIn = self.numZIn / self.numXIn
            self.yIn = y[0:self.numZIn:self.numXIn]
            self.zIn = numpy.reshape(z, (self.numYIn, self.numXIn))
        else:
            # Ordered by columns.
            self.numYIn = max(numpy.argmax(y) + 1, numpy.argmin(y) + 1)
            self.yIn = y[0:self.numYIn]
            self.numXIn = self.numZIn / self.numYIn
            self.xIn = x[0:self.numZIn:self.numYIn]
            self.ZIn = numpy.transpose(numpy.reshape(z, (self.numXIn, self.numYIn)))

        if (self.xIn[0] > self.xIn[1]):
            self.xIn = numpy.flipud(self.xIn)
            self.zIn = numpy.fliplr(self.zIn)
        if (self.yIn[0] > self.yIn[1]):
            self.yIn = numpy.flipud(self.yIn)
            self.zIn = numpy.flipud(self.zIn)

        return

    def _resampleDem(self):
        """
        Resample DEM in regions outside specified region.
        """

        # Find indices corresponding to min and max values.
        xMinInd = numpy.searchsorted(self.xIn, self.xMin, side='left')
        xMaxInd = numpy.searchsorted(self.xIn, self.xMax, side='right')
        yMinInd = numpy.searchsorted(self.yIn, self.yMin, side='left')
        yMaxInd = numpy.searchsorted(self.yIn, self.yMax, side='right')

        # Determine indices to extract.
        xIndex = xMinInd
        xIndices = []
        count = 0
        while (xIndex > 0):
            count += 1
            xIndex -= count * self.skipInterval
            xIndices.append(xIndex)
        numIndices = len(xIndices)
        if (numIndices != 0):
            if (xIndices[numIndices - 1] != 0):
                if (xIndices[numIndices - 1] < 0):
                    xIndices[numIndices - 1] = 0
                elif (xIndices[numIndices - 1] > 0):
                    xIndices.append(0)
        xIndices.reverse()
        xIndices.extend(range(xMinInd, min(xMaxInd + 1, self.numXIn)))

        xIndex = xMaxInd
        count = 0
        while (xIndex < self.numXIn - 1):
            count += 1
            xIndex += count * self.skipInterval
            xIndices.append(xIndex)
        numIndices = len(xIndices)
        if (xIndices[numIndices - 1] != self.numXIn - 1):
            if (xIndices[numIndices - 1] > self.numXIn - 1):
                xIndices[numIndices - 1] = self.numXIn - 1
            elif (xIndices[numIndices - 1] < self.numXIn - 1):
                xIndices.append(self.numXIn - 1)

        yIndex = yMinInd
        yIndices = []
        count = 0
        while (yIndex > 0):
            count += 1
            yIndex -= count * self.skipInterval
            yIndices.append(yIndex)
        numIndices = len(yIndices)
        if (numIndices != 0):
            if (yIndices[numIndices - 1] != 0):
                if (yIndices[numIndices - 1] < 0):
                    yIndices[numIndices - 1] = 0
                elif (yIndices[numIndices - 1] > 0):
                    yIndices.append(0)
        yIndices.reverse()
        yIndices.extend(range(yMinInd, min(yMaxInd + 1, self.numYIn)))

        yIndex = yMaxInd
        count = 0
        while (yIndex < self.numYIn - 1):
            count += 1
            yIndex += count * self.skipInterval
            yIndices.append(yIndex)
        numIndices = len(yIndices)
        if (yIndices[numIndices - 1] != self.numYIn - 1):
            if (yIndices[numIndices - 1] > self.numYIn - 1):
                yIndices[numIndices - 1] = self.numYIn - 1
            elif (yIndices[numIndices - 1] < self.numYIn - 1):
                yIndices.append(self.numYIn - 1)

        # Get output values.
        self.xOut = self.xIn[xIndices]
        self.yOut = self.yIn[yIndices]
        zTmp = numpy.take(self.zIn, yIndices, axis=0)
        self.zOut = numpy.take(zTmp, xIndices, axis=1)

        # Compute output sizes.
        self.numXOut = len(self.xOut)
        self.numYOut = len(self.yOut)
        self.numZOut = self.numXOut * self.numYOut

        # Smooth z-values that lie outside specified range.
        for xIndex in range(self.numXOut):
            xIndexIn = xIndices[xIndex]
            xCoord = self.xOut[xIndex]
            if (xIndex == 0):
                xDiff = xIndices[xIndex + 1] - xIndexIn
            elif (xIndex == self.numXOut - 1):
                xDiff = xIndexIn - xIndices[xIndex - 1]
            else:
                xDiff = max((xIndices[xIndex + 1] - xIndexIn),
                            (xIndexIn - xIndices[xIndex - 1]))
            for yIndex in range(self.numYOut):
                yIndexIn = yIndices[yIndex]
                yCoord = self.yOut[yIndex]
                if (yIndex == 0):
                    yDiff = yIndices[yIndex + 1] - yIndexIn
                elif (yIndex == self.numYOut - 1):
                    yDiff = yIndexIn - yIndices[yIndex - 1]
                else:
                    yDiff = max((yIndices[yIndex + 1] - yIndexIn),
                                (yIndexIn - yIndices[yIndex - 1]))
                if (xCoord < self.xMin or xCoord > self.xMax or
                        yCoord < self.yMin or yCoord > self.yMax):
                    self.zOut[yIndex, xIndex] = self._zavg(xDiff, yDiff,
                                                           xIndexIn, yIndexIn)

        return

    def _zavg(self, xDiff, yDiff, xIndexIn, yIndexIn):
        """
        Spatially average elevation values.
        """
        maxDiff = max(xDiff, yDiff)
        windowSize = maxDiff // 2
        xMin = max(xIndexIn - windowSize, 0)
        xMax = min(xIndexIn + windowSize + 1, self.numXIn - 1)
        yMin = max(yIndexIn - windowSize, 0)
        yMax = min(yIndexIn + windowSize + 1, self.numYIn - 1)
        zWindow = self.zIn[yMin:yMax, xMin:xMax]
        xTmp = numpy.fabs(numpy.arange(xMin - xIndexIn, xMax - xIndexIn))
        xWeight = numpy.fabs(xTmp - numpy.amax(xTmp)) + 1
        yTmp = numpy.fabs(numpy.arange(yMin - yIndexIn, yMax - yIndexIn))
        yWeight = numpy.fabs(yTmp - numpy.amax(yTmp)) + 1

        weight = numpy.outer(yWeight, xWeight)

        zAvg = numpy.average(zWindow, weights=weight)
        return zAvg

    def _writeCubitJournals(self):
        """
        Writes lines of the DEM as sets of Cubit journal files and write master
        journal file.
        """

        numWidth = 4
        fmt = " location %15.11e %15.11e %15.11e"
        newLine = "\n"
        masterPref = "playback '"
        separator = "# ----------------------------------------------------------\n"
        comment = "#"

        # Write master journal file.
        playCmd = \
            "reset\n" + \
            "# This is a simple Cubit journal file to create an ACIS\n" + \
            "# NURBS surface from a set of intersecting lines.\n" + \
            comment + newLine + separator + \
            "# Create u-lines and v-lines, then create a net surface.\n" + \
            separator + \
            "playback " + "'" + self.uLineJournal + "'" + "\n" + \
            "playback " + "'" + self.vLineJournal + "'" + "\n"
        m = open(self.masterJournal, 'w')
        m.write(playCmd)
        c1 = 1
        c2 = self.numYOut
        c3 = c2 + 1
        c4 = c2 + self.numXOut
        netCmd = "create surface net u curve " + repr(c1) + " to " + repr(c2) + \
                 " v curve " + repr(c3) + " to " + repr(c4) + newLine
        m.write(netCmd)
        comment2 = comment + newLine + separator + \
            "# Delete curves and any extra vertices.\n" + \
            separator
        m.write(comment2)
        delCmd = "delete curve " + repr(c1) + " to " + repr(c4) + newLine + \
                 "delete vertex all\n"
        m.write(delCmd)
        comment3 = comment + newLine + separator + \
            "# Export binary ACIS file.\n" + separator
        m.write(comment3)
        exportCmd = "export Acis '" + self.acisFilename + "'\n"
        m.write(exportCmd)
        m.close()

        # Write out u (east) lines.
        um = open(self.uLineJournal, 'w')
        for row in range(self.numYOut):
            y = self.yOut[row]
            uString = repr(row + 1).rjust(numWidth, '0')
            outputFileName = self.uLinePrefix + "_u" + uString + ".jou"
            masterString = masterPref + outputFileName + "'" + newLine
            um.write(masterString)
            u = open(outputFileName, 'w')
            for column in range(self.numXOut):
                point = (self.xOut[column], y, self.zOut[row, column])
                u.write("create vertex x %10.5e y %10.5e z %10.5e\n" %
                        (point[0], point[1], point[2]))
                if 0 == column:
                    u.write("${idBeg=Id('vertex')}\n")
            u.write("${idEnd=Id('vertex')}\n")
            u.write("create curve spline vertex {idBeg} to {idEnd} delete\n")

            u.close()

        um.close()
        print("Number of u-lines = " + repr(self.numYOut))

        # Write out v (north) lines.
        vm = open(self.vLineJournal, 'w')
        for column in range(self.numXOut):
            x = self.xOut[column]
            vString = repr(column + 1).rjust(numWidth, '0')
            outputFileName = self.vLinePrefix + "_v" + vString + ".jou"
            masterString = masterPref + outputFileName + "'" + newLine
            vm.write(masterString)
            v = open(outputFileName, 'w')
            for row in range(self.numYOut):
                point = (x, self.yOut[row], self.zOut[row, column])
                v.write("create vertex x %10.5e y %10.5e z %10.5e\n" %
                        (point[0], point[1], point[2]))
                if 0 == row:
                    v.write("${idBeg=Id('vertex')}\n")
            v.write("${idEnd=Id('vertex')}\n")
            v.write("create curve spline vertex {idBeg} to {idEnd} delete\n")
            v.close()

        vm.close()
        print("Number of v-lines = " + repr(self.numXOut))

        return

    def _writeDemVtk(self):
        """
        Write DEM as a rectilinear grid VTK file with z-values as point data.
        """
        zDim = 1
        v = open(self.vtkOutputFile, 'w')
        v.write('# vtk DataFile Version 2.0\n')
        v.write('Resampled DEM\n')
        v.write('ASCII\n')
        v.write('DATASET RECTILINEAR_GRID\n')
        dimString = 'DIMENSIONS ' + str(self.numXOut) + ' ' + str(self.numYOut) + \
                    ' ' + str(zDim) + '\n'
        v.write(dimString)

        xString = 'X_COORDINATES ' + str(self.numXOut) + ' double\n'
        v.write(xString)
        for point in range(self.numXOut):
            v.write("%15.11e  " % self.xOut[point])
            if ((point + 1) % 5 == 0):
                v.write("\n")

        yString = '\nY_COORDINATES ' + str(self.numYOut) + ' double\n'
        v.write(yString)
        for point in range(self.numYOut):
            v.write("%15.11e  " % self.yOut[point])
            if ((point + 1) % 5 == 0):
                v.write("\n")

        zString = '\nZ_COORDINATES ' + str(zDim) + ' double\n'
        v.write(zString)
        v.write('0.0\n')

        zString1 = 'POINT_DATA ' + str(self.numZOut) + '\n'
        v.write(zString1)
        zString2 = 'SCALARS elevation double 1\n'
        v.write(zString2)
        zString3 = 'LOOKUP_TABLE default\n'
        v.write(zString3)
        for yPoint in range(self.numYOut):
            for xPoint in range(self.numXOut):
                v.write("%15.11e  " % self.zOut[yPoint, xPoint])
                if ((xPoint + 1) % 5 == 0):
                    v.write("\n")

        v.close()
        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = Dem2Lines()
    app.run()

# End of file
