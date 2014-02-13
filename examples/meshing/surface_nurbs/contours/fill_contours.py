#!/usr/bin/env python

## @file fill_contours.py

## @brief Python application to read a set of partial contours lying between
## two complete contours, and fill in missing values. It is assumed that we
## want approximately equal spacing in the y-direction.

import math
import numpy
import pdb

from pyre.applications.Script import Script as Application

class FillContours(Application):
  """
  Python application to read a set of partial contours lying between
  two complete contours, and fill in missing values. It is assumed that we
  want approximately equal spacing in the y-direction.
  """

  class Inventory(Application.Inventory):
    """
    Python object for managing FillContours facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FillContours facilities and properties.
    ##
    ## \b Properties
    ## @li \b in_file Input file containing all contours.
    ## @li \b out_file Output file containing altered contours.
    ## @li \b interval_size Desired y-interval.
    ## @li \b min_check Minimum y-value to check for interpolation gaps.
    ## @li \b max_check Maximum y-value to check for interpolation gaps.

    import pyre.inventory
    from pyre.units.angle import degree
    
    inFile = pyre.inventory.str("in_file", default="contours_in.txt")
    inFile.meta['tip'] = "Input file containing original contours."
    
    outFile = pyre.inventory.str("out_file", default="contours_out.txt")
    outFile.meta['tip'] = "Output file containing altered contours."
    
    intervalSize = pyre.inventory.float("interval_size", default=1.0)
    intervalSize.meta['tip'] = "Desired y-interval size."
    
    minCheck = pyre.inventory.float("min_check", default=-1.0e8)
    minCheck.meta['tip'] = "Minimum y-value to check for interpolation gaps."
    
    maxCheck = pyre.inventory.float("max_check", default=1.0e8)
    maxCheck.meta['tip'] = "Maximum y-value to check for interpolation gaps."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fill_contours"):
    Application.__init__(self, name)
    self.numContours = 0
    self.contourValues = []
    self.numPointsPerContour = []
    self.contoursIn = []
    self.contoursOut = []
    self.spaceDim = 3
    self.minYValue = 0.0
    self.maxYValue = 0.0

    return


  def main(self):
    # pdb.set_trace()
    self._readContours()
    self._fillContours()
    self._writeContours()

    return
  

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)

    # Filenames
    self.inFile = self.inventory.inFile
    self.outFile = self.inventory.outFile

    # Parameters
    self.intervalSize = self.inventory.intervalSize
    self.minCheck = self.inventory.minCheck
    self.maxCheck = self.inventory.maxCheck

    return


  def _readContours(self):
    """
    Read points defining contours and break them into separate arrays.
    """
    points = []
    f = open(self.inFile, 'r')
    totalPoints = 0
    contourPrev = -9999.0
    numPointsInContour = 0

    for line in f:
      data = line.split()
      point = [float(data[0]), float(data[1]), float(data[2])]
      points.append(point)
      totalPoints += 1
      numPointsInContour += 1
      if (totalPoints == 1):
        self.numContours += 1
        self.contourValues.append(float(data[2]))
      elif (float(data[2]) != self.contourValues[self.numContours - 1]):
        self.numPointsPerContour.append(numPointsInContour - 1)
        self.contourValues.append(float(data[2]))
        self.numContours += 1
        numPointsInContour = 1

    self.numPointsPerContour.append(numPointsInContour)

    f.close()

    contourIndex = 0
    
    for contour in range(self.numContours):
      self.contoursIn.append(
        numpy.array(
        points[contourIndex:contourIndex+self.numPointsPerContour[contour]],
        dtype=numpy.float64).reshape(self.numPointsPerContour[contour],
                                    self.spaceDim))
      contourIndex += self.numPointsPerContour[contour]
      yContour = self.contoursIn[contour][:,1]
      contourMax = yContour.max()
      contourMin = yContour.min()
      self.maxYValue = max(self.maxYValue, contourMax)
      self.minYValue = min(self.minYValue, contourMin)

    return

  
  def _fillContours(self):
    """
    Add points to contours with missing pieces.
    """
    # We assume that the first and last contours are OK.
    indMinContour = 0
    cMinContour = self.contourValues[indMinContour]
    xMinContour = self.contoursIn[indMinContour][:,0]
    yMinContour = self.contoursIn[indMinContour][:,1]
    indMaxContour = self.numContours - 1
    cMaxContour = self.contourValues[indMaxContour]
    xMaxContour = self.contoursIn[indMaxContour][:,0]
    yMaxContour = self.contoursIn[indMaxContour][:,1]

    self.contoursOut = [self.contoursIn[indMinContour]]
    
    for contour in range(1, self.numContours - 1):
      cContour = self.contourValues[contour]
      xContour = self.contoursIn[contour][:,0]
      yContour = self.contoursIn[contour][:,1]
      contourYMin = yContour.min()
      contourYMax = yContour.max()
      point = numpy.array(([0.0, self.minYValue, cContour]),
                          dtype=numpy.float64, ndmin=2)

      # If contour doesn't contain minimum y-value, add it.
      startSearch = 0
      if (contourYMin == self.minYValue):
        point[0,0] = xContour[0]
        startSearch = 1
      else:
        point[0,0] = self._findValue(contour, self.minYValue)
      contoursOut = point
      yPrev = self.minYValue

      # Look for gaps larger than the required interval, and fill in as needed.
      for pointIndex in range(startSearch, self.numPointsPerContour[contour]):
        yVal = yContour[pointIndex]
        yDiff = yVal - yPrev
        if(yDiff > self.intervalSize):
          numNewContours = int(math.ceil(yDiff/self.intervalSize))
          contInt = yDiff/numNewContours
          newPoint = numpy.empty((1,3), dtype=numpy.float64)
          newPoint[0,2] = cContour
          for newContour in range(numNewContours - 1):
            yNew = yPrev + (newContour + 1) * contInt
            newPoint[0,1] = yNew
            newPoint[0,0] = self._findValue(contour, yNew)
            contoursOut = numpy.append(contoursOut, newPoint, axis=0)
        
        existingPoint = numpy.array(([xContour[pointIndex],
                                      yContour[pointIndex],
                                      cContour]),
                                    dtype=numpy.float64, ndmin=2)
        contoursOut = numpy.append(contoursOut, existingPoint, axis=0)
        yPrev = yVal
        
      # Look for gap above maximum.
      if (contourYMax != self.maxYValue):
        yDiff = self.maxYValue - contourYMax
        numNewContours = int(math.ceil(yDiff/self.intervalSize))
        contInt = yDiff/numNewContours
        newPoint = numpy.empty((1,3), dtype=numpy.float64)
        newPoint[0,2] = cContour
        for newContour in range(numNewContours):
          yNew = contourYMax + (newContour + 1) * contInt
          newPoint[0,1] = yNew
          newPoint[0,0] = self._findValue(contour, yNew)
          contoursOut = numpy.append(contoursOut, newPoint, axis=0)

      self.contoursOut.append(contoursOut) 

    # Append last contour
    self.contoursOut.append(self.contoursIn[indMaxContour])

    return

      
  def _findValue(self, contour, yVal):
    """
    Search through existing contours and interpolate to get x-coordinate of
    requested point.
    """
    # pdb.set_trace()
    # We use the next lowest contour as our lower value, since these have
    # already been created.
    contourLow = contour - 1
    xLows = self.contoursOut[contourLow][:,0]
    yLows = self.contoursOut[contourLow][:,1]
    xLow = numpy.interp(yVal, yLows, xLows)
    contourLowValue = self.contourValues[contourLow]

    # Find the first contour above this one that doesn't have too large a gap.
    xHigh = 0.0
    contourHigh = 0
    for testContour in range(contour + 1, self.numContours):
      yHighs = self.contoursIn[testContour][:,1]
      numContourPoints = self.numPointsPerContour[testContour]
      testIndex = 0
      if (yVal == self.minYValue):
        testIndex = numpy.searchsorted(yHighs, yVal, side='right')
      else:
        testIndex = numpy.searchsorted(yHighs, yVal, side='left')
      if (testIndex > 0 and testIndex < numContourPoints):
        yBelow = yHighs[testIndex - 1]
        yAbove = yHighs[testIndex]
        yDiff = yAbove - yBelow
        xHighs = self.contoursIn[testContour][:,0]
        if (yDiff < self.intervalSize or \
            yVal < self.minCheck or yVal > self.maxCheck):
          xHighs = self.contoursIn[testContour] [:,0]
          xHigh = numpy.interp(yVal, yHighs, xHighs)
          contourHigh = testContour
          break
#        if (yVal < self.minCheck):
#          xHigh = xHighs[0]
#          contourHigh = testContour
#          break
#        elif (yVal > self.maxCheck):
#          xHigh = xHighs[numContourPoints - 1]
#          contourHigh = testContour
#          break
#        elif (yDiff < self.intervalSize):
#          xHigh = numpy.interp(yVal, yHighs, xHighs)
#          contourHigh = testContour
#          break

    # Given the values above and below, interpolate to find the correct x-value.
    contourHighValue = self.contourValues[contourHigh]
    contourValue = self.contourValues[contour]
    cDiffHighLow = contourHighValue - contourLowValue
    xDiffHighLow = xHigh - xLow
    ratioXC = xDiffHighLow/cDiffHighLow
    xInterp = xLow + (contourValue - contourLowValue) * ratioXC

    return xInterp


  def _writeContours(self):
    """
    Write out the revised contours to a single file.
    """
    output = self.contoursOut[0]
    for contour in range(1, self.numContours):
      output = numpy.append(output, self.contoursOut[contour], axis=0)

    numpy.savetxt(self.outFile, output)

    return

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = FillContours()
  app.run()

# End of file

