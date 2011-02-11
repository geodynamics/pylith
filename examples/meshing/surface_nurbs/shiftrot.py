#!/usr/bin/env python

## @file shiftrot.py

## @brief Python application to take a set of (x,y,z) points from an ASCII file,
## shift them in the (x,y) plane by a given amount, and then rotate them about
## the origin by a given amount in the (x,y) plane (shifting occurs in the old
## coordinate system).
## It is assumed that the input coordinates and the shift values have the same
## units.

import math
import pdb

from pyre.applications.Script import Script as Application

class ShiftRot(Application):
  """
  Python application to take a set of (x,y,z) points from an ASCII file,
  shift them in the (x,y) plane by a given amount, and then rotate them about
  the origin by a given amount in the (x,y) plane (shifting occurs in the old
  coordinate system).
  It is assumed that the input coordinates and the shift values have the same
  units.
  """

  class Inventory(Application.Inventory):
    """
    Python object for managing ShiftRot facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing ShiftRot facilities and properties.
    ##
    ## \b Properties
    ## @li \b in_file Input file containing original coordinates.
    ## @li \b out_file Output file containing shifted and rotated coordinates.
    ## @li \b x_shift Amount by which to shift coordinates in the x-direction.
    ## @li \b y_shift Amount by which to shift coordinates in the y-direction.
    ## @li \b rot_angle Amount by which to rotate coordinates.
    ## @li \b scale_factor Amount by which to scale coordinates.

    import pyre.inventory
    from pyre.units.angle import degree
    
    inFile = pyre.inventory.str("in_file", default="in_file.txt")
    inFile.meta['tip'] = "Input file containing original coordinates."
    
    outFile = pyre.inventory.str("out_file", default="out_file.txt")
    outFile.meta['tip'] = "Output file with shifted and rotated coordoutates."
    
    xShift = pyre.inventory.float("x_shift", default=0.0)
    xShift.meta['tip'] = "X-shift to apply to coordinates."
    
    yShift = pyre.inventory.float("y_shift", default=0.0)
    yShift.meta['tip'] = "Y-shift to apply to coordinates."
    
    rotAngle = pyre.inventory.dimensional("rot_angle", default=0.0*degree)
    rotAngle.meta['tip'] = "Amount by which to rotate coordinates."
    
    scaleFactor = pyre.inventory.float("scale_factor", default=1.0)
    scaleFactor.meta['tip'] = "Scaling factor to apply to results."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="shiftrot"):
    Application.__init__(self, name)

    return


  def main(self):
    # pdb.set_trace()
    self._transformCoords()

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
    self.xShift = self.inventory.xShift
    self.yShift = self.inventory.yShift
    self.rotAngle = self.inventory.rotAngle.value
    self.scaleFactor = self.inventory.scaleFactor

    return


  def _transformCoords(self):
    """
    Rotate, shift, and scale coordinates.
    """
    # pdb.set_trace()
    cosRot = math.cos(self.rotAngle)
    sinRot = math.sin(self.rotAngle)

    f = open(self.inFile, 'r')
    g = open(self.outFile, 'w')

    for line in f:
      data = line.split()
      x = float(data[0])
      y = float(data[1])
      z = float(data[2])
      xTrans = self.scaleFactor * ( (x + self.xShift) * cosRot + \
                                    (y + self.yShift) * sinRot)
      yTrans = self.scaleFactor * (-(x + self.xShift) * sinRot + \
                                   (y + self.yShift) * cosRot)
      g.write('%.12f' % xTrans)
      g.write('  %.12f' % yTrans)
      g.write('  %.12f' % z)
      g.write( '\n')

    f.close()
    g.close()

    return

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = ShiftRot()
  app.run()

# End of file

