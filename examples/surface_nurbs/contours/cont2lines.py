#!/usr/bin/env python

## @file cont2lines.py

## @brief Python application to read a set of contours from a file
## and output them as lines to be used with a Cubit script.
## A control script is also written than will make the lines into a
## skin surface.

import pdb

from pyre.applications.Script import Script as Application

class Cont2Lines(Application):
  """
  Python application to read a set of contours from a file
  and output them as lines to be used with a Cubit script.
  A control script is also written than will make the lines into a
  skin surface.
  """

  class Inventory(Application.Inventory):
    """
    Python object for managing Cont2Lines facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Cont2Lines facilities and properties.
    ##
    ## \b Properties
    ## @li \b in_file Input file containing all contours.
    ## @li \b out_root Root filename for each journal file defining a line.
    ## @li \b journal_out Name of controlling journal file.
    ## @li \b surface_out Name of output surface (.sab file).

    import pyre.inventory
    
    inFile = pyre.inventory.str("in_file", default="contours_in.txt")
    inFile.meta['tip'] = "Input file containing original contours."
    
    outRoot = pyre.inventory.str("out_root", default="contours_out")
    outRoot.meta['tip'] = "Root filename for each line definition journal file."
    
    journalOut = pyre.inventory.str("journal_out", default="mksurf.jou")
    journalOut.meta['tip'] = "Name of controlling journal file."
    
    surfaceOut = pyre.inventory.str("surface_out", default="surface.sab")
    surfaceOut.meta['tip'] = "Name of output surface (.sab file)."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="cont2lines"):
    Application.__init__(self, name)

    return


  def main(self):
    # pdb.set_trace()
    self._readContours()

    return
  

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)

    # Filenames
    self.inFile = self.inventory.inFile
    self.outRoot = self.inventory.outRoot
    self.journalOut = self.inventory.journalOut
    self.surfaceOut = self.inventory.surfaceOut

    return


  def _readContours(self):
    """
    Read points defining contours and output each one as a journal file.
    """
    # Read contours and determine how many points in each.
    f = open(self.inFile, 'r')
    totalPoints = 0
    contourPrev = -9999.0
    contourIndex = 0
    fmt = " location %15.11e %15.11e %15.11e"
    indexWidth = 4
    fileOut = self.outRoot + '_c' + \
              repr(contourIndex).rjust(indexWidth, '0') + '.jou'
    c = open(fileOut, 'w')
    j = open(self.journalOut, 'w')
    j.write('reset\n')
    j.write('# --------------------------------------------------------\n')
    j.write('# Create a spline for each contour line.\n')
    j.write('# --------------------------------------------------------\n')


    for line in f:
      data = line.split()
      point = [float(data[0]), float(data[1]), float(data[2])]
      totalPoints += 1
      if (totalPoints == 1):
        c.write('create curve spline')
        c.write(fmt % (point[0], point[1], point[2]))
        j.write('playback \'' + fileOut + '\'\n')
      elif (float(data[2]) != contourPrev):
        contourIndex += 1
        c.write('\n')
        c.close()
        fileOut = self.outRoot + '_c' + \
                  repr(contourIndex).rjust(indexWidth, '0') + '.jou'
        c = open(fileOut, 'w')
        c.write('create curve spline')
        c.write(fmt % (point[0], point[1], point[2]))
        j.write('playback \'' + fileOut + '\'\n')
      else:
        c.write(fmt % (point[0], point[1], point[2]))
      contourPrev = float(data[2])

    j.write('# --------------------------------------------------------\n')
    j.write('# Create skin surface from contours and then delete contours.\n')
    j.write('# --------------------------------------------------------\n')
    j.write('create surface skin curve all\n')
    firstCurve = 1
    lastCurve = contourIndex + 1
    j.write('delete curve %d to %d\n' % (firstCurve, lastCurve))
    j.write('# --------------------------------------------------------\n')
    j.write('# Export ACIS file.\n')
    j.write('# --------------------------------------------------------\n')
    j.write('export Acis \'' + self.surfaceOut + '\'\n')
    c.close()
    j.close()

    return

  
# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = Cont2Lines()
  app.run()

# End of file

