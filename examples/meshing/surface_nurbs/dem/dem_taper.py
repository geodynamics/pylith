#!/usr/bin/env python

## @file dem_taper.py

## @brief Python application to read an ASCII (x,y,z) representation of a DEM
## and create a set of lines suitable for use in Cubit. The original DEM is
## tapered at the edges to a constant value.

import math
import numpy
import pdb

from pyre.applications.Script import Script as Application

class DemTaper(Application):
  """
  Python application to read an ASCII (x,y,z) representation and
  create a set of lines suitable for use in Cubit. The original DEM is
  tapered at the edges to a constant value.
  The DEM is assumed to be ordered by rows (left to right).
  """

  class Inventory(Application.Inventory):
    """
    Python object for managing DemTaper facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing DemTaper facilities and properties.
    ##
    ## \b Properties
    ## @li \b input_dem Input DEM file (ASCII).
    ## @li \b vtk_output_file Filename of VTK output file.
    ## @li \b u_line_prefix Prefix for u (east) output lines.
    ## @li \b u_line_journal Output journal file for u_lines.
    ## @li \b v_line_prefix Prefix for v (north) output lines.
    ## @li \b v_line_journal Output journal file for v_lines.
    ## @li \b master_journal Output master journal file.
    ## @li \b acis_filename Name of ACIS output file created by Cubit.
    ## @li \b taper_len Number of pixels over which to taper DEM.
    ## @li \b edge_value_source How to determine DEM values on outer edges.
    ## @li \b edge_value User-specified edge value.
    ## @li \b interpolation_type Interpolation method to use.

    import pyre.inventory
    from pyre.units.angle import degree
    
    inputDem = pyre.inventory.str("input_dem", default="DEM.txt")
    inputDem.meta['tip'] = "Input DEM file."
    
    vtkOutputFile = pyre.inventory.str("vtk_output_file", default="DEM.vtk")
    vtkOutputFile.meta['tip'] = "VTK output file."
    
    uLinePrefix = pyre.inventory.str("u_line_prefix", default="u_lines")
    uLinePrefix.meta['tip'] = "Prefix for u (east) lines."
    
    uLineJournal = pyre.inventory.str("u_line_journal", default="u_lines.jou")
    uLineJournal.meta['tip'] = "Output journal file for u (east) lines."
    
    vLinePrefix = pyre.inventory.str("v_line_prefix", default="v_lines")
    vLinePrefix.meta['tip'] = "Prefix for v (north) lines."
    
    vLineJournal = pyre.inventory.str("v_line_journal", default="v_lines.jou")
    vLineJournal.meta['tip'] = "Output journal file for v (north) lines."
    
    masterJournal = pyre.inventory.str("master_journal", default="master.jou")
    masterJournal.meta['tip'] = "Output master journal file."
    
    acisFilename = pyre.inventory.str("acis_filename", default="topo.sab")
    acisFilename.meta['tip'] = "Name of ACIS output file created by Cubit."
    
    taperLen = pyre.inventory.int("taper_len", default=10)
    taperLen.meta['tip'] = "Number of pixels over which to taper DEM."
    
    edgeValueSource = pyre.inventory.str("edge_value_source",
      default="mean",
      validator=pyre.inventory.choice(["mean","user_specified"]))
    edgeValueSource.meta['tip'] = "How to determine DEM values on outer edges."
    
    edgeValue = pyre.inventory.int("edge_value", default=0)
    edgeValue.meta['tip'] = "User-specified value for DEM edges."
    
    interpolationType = pyre.inventory.str("interpolation_type",
      default="linear",
      validator=pyre.inventory.choice(["nearest","linear","cubic"]))
    interpolationType.meta['tip'] = "Interpolation method to use."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="dem_taper"):
    Application.__init__(self, name)
    self.xDim = 0
    self.yDim = 0
    self.zDim = 0
    self.meanZ = 0
    self.keepIndicesX = []
    self.replaceIndicesX = []
    self.keepIndicesY = []
    self.replaceIndicesY = []
    self.xIn = None
    self.yIn = None
    self.zIn = None
    self.xLine = None
    self.yLine = None
    self.valuesOut = None

    return


  def main(self):
    # pdb.set_trace()
    self._readDem()
    self._taperDem()
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
    self.uLinePrefix = self.inventory.uLinePrefix
    self.uLineJournal = self.inventory.uLineJournal
    self.vLinePrefix = self.inventory.vLinePrefix
    self.vLineJournal = self.inventory.vLineJournal
    self.masterJournal = self.inventory.masterJournal
    self.acisFilename = self.inventory.acisFilename

    # Parameters
    self.taperLen = self.inventory.taperLen
    self.edgeValueSource = self.inventory.edgeValueSource
    self.edgeValue = self.inventory.edgeValue
    self.interpolationType = self.inventory.interpolationType

    return


  def _readDem(self):
    """
    Read coordinates defining DEM and create vectors of x, y, and z values.
    """

    # Load each coordinate as a numpy array.
    x, y, z = numpy.loadtxt(self.inputDem, dtype=numpy.float64, unpack=True)

    self.zDim = len(z)
    if (y[0] == y[1]):
      # Ordered by rows.
      self.xDim = max(numpy.argmax(x) + 1, numpy.argmin(x) + 1)
      self.xLine = x[0:self.xDim]
      self.yDim = self.zDim/self.xDim
      self.yLine = y[0:self.zDim:self.xDim]
      self.xIn = numpy.reshape(x, (self.yDim, self.xDim))
      self.yIn = numpy.reshape(y, (self.yDim, self.xDim))
      self.zIn = numpy.reshape(z, (self.yDim, self.xDim))
    else:
      # Ordered by columns.
      self.yDim = max(numpy.argmax(y) + 1, numpy.argmin(y) + 1)
      self.yLine = y[0:self.yDim]
      self.xDim = self.zDim/self.yDim
      self.xLine = x[0:self.zDim:self.yDim]
      self.xIn = numpy.transpose(numpy.reshape(x, (self.xDim, self.yDim)))
      self.yIn = numpy.transpose(numpy.reshape(y, (self.xDim, self.yDim)))
      self.zIn = numpy.transpose(numpy.reshape(z, (self.xDim, self.yDim)))

    # Flip everything so that it is ordered left-to-right and bottom-to-top.
    if (self.xLine[0] > self.xLine[1]):
      self.xLine = numpy.flipud(self.xLine)
      self.xIn = numpy.fliplr(self.xIn)
      self.yIn = numpy.fliplr(self.yIn)
      self.zIn = numpy.fliplr(self.zIn)
    if (self.yLine[0] > self.yLine[1]):
      self.yLine = numpy.flipud(self.yLine)
      self.xIn = numpy.flipud(self.xIn)
      self.yIn = numpy.flipud(self.yIn)
      self.zIn = numpy.flipud(self.zIn)

    return

  
  def _taperDem(self):
    """
    Taper DEM over specified region.
    """

    import matplotlib.mlab

    # Determine regions over which to replace existing values.
    xLeftOuter = 0
    xLeftInner = xLeftOuter + self.taperLen
    xRightOuter = self.xDim - 1
    xRightInner = xRightOuter - self.taperLen
    yBottomOuter = 0
    yBottomInner = yBottomOuter + self.taperLen
    yTopOuter = self.yDim - 1
    yTopInner = yTopOuter - self.taperLen
    self.keepIndicesX = range(xLeftInner + 1, xRightInner - 1)
    self.replaceIndicesX = range(xLeftOuter, xLeftInner) + \
                           range(xRightInner, xRightOuter)
    self.keepIndicesY = range(yBottomInner + 1, yTopInner - 1)
    self.replaceIndicesY = range(yBottomOuter, yBottomInner) + \
                           range(yTopInner, yTopOuter)

    # Get edge values.
    xEdge1 = self.xIn[yBottomOuter, :]
    yEdge1 = self.yIn[yBottomOuter, :]
    zEdge1 = self.zIn[yBottomOuter, :]
    xEdge2 = self.xIn[yTopOuter, :]
    yEdge2 = self.yIn[yTopOuter, :]
    zEdge2 = self.zIn[yTopOuter, :]
    xEdge3 = self.xIn[yBottomOuter+1:yTopOuter, xLeftOuter]
    yEdge3 = self.yIn[yBottomOuter+1:yTopOuter, xLeftOuter]
    zEdge3 = self.zIn[yBottomOuter+1:yTopOuter, xLeftOuter]
    xEdge4 = self.xIn[yBottomOuter+1:yTopOuter, xRightOuter]
    yEdge4 = self.yIn[yBottomOuter+1:yTopOuter, xRightOuter]
    zEdge4 = self.zIn[yBottomOuter+1:yTopOuter, xRightOuter]
    zEdge = numpy.concatenate((zEdge1, zEdge2, zEdge3, zEdge4))
    self.zMean = int(numpy.mean(zEdge))
    edgeValue = self.zMean
    if (self.edgeValueSource == "user_specified"):
      edgeValue = self.edgeValue

    print "Mean elevation around DEM edges:  %i" % self.zMean
    
    # Create arrays with taper values missing.
    xTmp = numpy.take(self.xIn, self.keepIndicesY, axis=0)
    xInner = numpy.take(xTmp, self.keepIndicesX, axis=1)
    xOuter = numpy.concatenate((xEdge1, xEdge2, xEdge3, xEdge4))
    xReduced = numpy.concatenate((xInner.flatten(), xOuter))
    yTmp = numpy.take(self.yIn, self.keepIndicesY, axis=0)
    yInner = numpy.take(yTmp, self.keepIndicesX, axis=1)
    yOuter = numpy.concatenate((yEdge1, yEdge2, yEdge3, yEdge4))
    yReduced = numpy.concatenate((yInner.flatten(), yOuter))
    zTmp = numpy.take(self.zIn, self.keepIndicesY, axis=0)
    zInner = numpy.take(zTmp, self.keepIndicesX, axis=1)
    zOuter = float(edgeValue) * numpy.ones_like(yOuter)
    zReduced = numpy.concatenate((zInner.flatten(), zOuter))

    # Interpolate missing values using selected method.
    self.valuesOut = matplotlib.mlab.griddata(xReduced, yReduced, zReduced,
                                              self.xIn, self.yIn)

    return
    
    
  def _writeCubitJournals(self):
    """
    Writes Cubit journal files to create a NURBS surface representing the DEM.
    """

    numWidth = 4
    fmt = " location %15.11e %15.11e %15.11e"
    newLine = "\n"
    masterPref = "playback '"
    separator = "# ------------------------------------------------------------"


    # Write out u (east) lines.
    um = open(self.uLineJournal, 'w')
    for row in range(self.yDim):
      y = self.yLine[row]
      uString = repr(row + 1).rjust(numWidth, '0')
      outputFileName = self.uLinePrefix + "_u" + uString + ".jou"
      masterString = masterPref + outputFileName + "'" + newLine
      um.write(masterString)
      u = open(outputFileName, 'w')
      u.write('create curve spline')
      for column in range(self.xDim):
        u.write(fmt % (self.xLine[column], y, self.valuesOut[row, column]))

      u.close()

    um.close()
    print "Number of u-lines = " + repr(self.yDim)

    # Write out v (north) lines.
    vm = open(self.vLineJournal, 'w')
    for column in range(self.xDim):
      x = self.xLine[column]
      vString = repr(column + 1).rjust(numWidth, '0')
      outputFileName = self.vLinePrefix + "_v" + vString + ".jou"
      masterString = masterPref + outputFileName + "'" + newLine
      vm.write(masterString)
      v = open(outputFileName, 'w')
      v.write('create curve spline')
      for row in range(self.yDim):
        v.write(fmt % (x, self.yLine[row], self.valuesOut[row, column]))

      v.close()

    vm.close()
    print "Number of v-lines = " + repr(self.xDim)

    # Write master journal file.
    m = open(self.masterJournal, 'w')
    m.write("reset" + newLine)
    m.write(separator + newLine)
    comment1 = "# Create u-lines and v-lines, and then create a net surface."
    m.write(comment1 + newLine)
    m.write(separator + newLine)
    m.write("playback '" + self.uLineJournal + "'" + newLine)
    m.write("playback '" + self.vLineJournal + "'" + newLine)

    uBegin = 1
    uEnd = self.yDim
    vBegin = uEnd + 1
    vEnd = vBegin + self.xDim - 1
    surfString = "create surface net u curve " + repr(uBegin) + " to " + \
                 repr(uEnd) + " v curve " + repr(vBegin) + " to " + \
                 repr(vEnd) + newLine
    m.write(surfString)
    m.write(newLine)

    m.write(separator + newLine)
    comment2 = "# Delete curves and any extra vertices." + newLine
    m.write(comment2)
    m.write(separator + newLine)
    curveDel = "delete curve " + repr(uBegin) + " to " + repr(vEnd) + newLine
    m.write(curveDel)
    m.write("delete vertex all" + newLine)
    m.write(newLine)
    
    m.write(separator + newLine)
    comment3 = "# Export binary ACIS file." + newLine
    m.write(comment3)
    m.write(separator + newLine)
    exportCmd = "export Acis '" + self.acisFilename + "'" + newLine
    m.write(exportCmd)
    m.close()

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
    dimString = 'DIMENSIONS ' + str(self.xDim) + ' ' + str(self.yDim) + \
                ' ' + str(zDim) + '\n'
    v.write(dimString)

    xString = 'X_COORDINATES ' + str(self.xDim) + ' double\n'
    v.write(xString)
    for point in range(self.xDim):
      v.write("%15.11e  " % self.xLine[point])
      if ((point + 1)%5 == 0):
        v.write("\n")

    yString = '\nY_COORDINATES ' + str(self.yDim) + ' double\n'
    v.write(yString)
    for point in range(self.yDim):
      v.write("%15.11e  " % self.yLine[point])
      if ((point + 1)%5 == 0):
        v.write("\n")

    zString = '\nZ_COORDINATES ' + str(zDim) + ' double\n'
    v.write(zString)
    v.write('0.0\n')

    zString1 = 'POINT_DATA ' + str(self.zDim) + '\n'
    v.write(zString1)
    zString2 = 'SCALARS elevation double 1\n'
    v.write(zString2)
    zString3 = 'LOOKUP_TABLE default\n'
    v.write(zString3)
    for yPoint in range(self.yDim):
      for xPoint in range(self.xDim):
        v.write("%15.11e  " % self.valuesOut[yPoint, xPoint])
        if ((xPoint + 1)%5 == 0):
          v.write("\n")

    v.close()      
    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = DemTaper()
  app.run()

# End of file

