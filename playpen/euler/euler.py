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
    ## Python object for managing Verifier facilities and properties.
    ##
    ## \b Properties
    ## @li \b points_file Filename of file containing point coordinates.
    ## @li \b points_datum Datum used for point coordinates.
    ## @li \b points_zone UTM zone used for point coordinates.
    ## @li \b points_offset_east East offset of point coords from UTM coords.
    ## @li \b points_offset_north North offset of point coords from UTM coords.
    ## @li \b euler_lat Latitude of Euler pole.
    ## @li \b euler_lon Longitude of Euler pole.
    ## @li \b euler_rot Rotation value for Euler pole (CCW positive).
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    pointsFile = pyre.inventory.str("points_file", default="points.def")
    pointsFile.meta['tip'] = "Filename of file containing point coordinates."

    pointsDatum = pyre.inventory.str("points_datum", default="WGS84")
    pointsDatum.meta['tip'] = "Datum used for point coordinates."

    pointsZone = pyre.inventory.str("points_zone", default="11")
    pointsZone.meta['tip'] = "UTM zone used for point coordinates."

    pointsOffsetEast = pyre.inventory.str("points_offset_east", default="0.0")
    pointsOffsetEast.meta['tip'] = "East offset of point coords from UTM coords."

    pointsOffsetNorth = pyre.inventory.str("points_offset_north", default="0.0")
    pointsOffsetNorth.meta['tip'] = "North offset of point coords from UTM coords."

    eulerLat = pyre.inventory.str("euler_lat", default="0.0")
    eulerLat.meta['tip'] = "Latitude of Euler pole."

    eulerLon = pyre.inventory.str("euler_lon", default="0.0")
    eulerLon.meta['tip'] = "Longitude of Euler pole."

    eulerRot = pyre.inventory.str("euler_rot", default="0.0")
    eulerRot.meta['tip'] = "Rotation of Euler pole (CCW positive)."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="euler"):
    Application.__init__(self, name)
    return


  def main(self):
    outputFiles = []
    for p in [1, 2]:
      basename = os.path.splitext(os.path.basename(self.cfgFile))[0] + \
                 '_p'+str(p)
      outputname = basename+'.vtk'
      outputFiles.append(basename+'_t0.vtk')
      cmd = '%s --nodes=%d %s --problem.formulation.output.output.filename=%s %s' % \
            (self.pylith, p, ' '.join(self.petscOptions), outputname,
             self.cfgFile)
      print 'Running %s' % cmd
      os.system(cmd)
    for file1, file2 in zip(outputFiles[:-1], outputFiles[1:]):
      self._compare(file1, file2)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.cfgFile = self.inventory.cfgFile
    self.pylith = self.inventory.pylith
    self.petscOptions = self.inventory.petscOptions
    return


  def _readVTK(self, filename):
    f = file(filename)
    found = False
    data = []
    for line in f.readlines():
      if self.exp.match(line):
        found = True
      if found and \
             not line.startswith('LOOKUP') and \
             not line.startswith('SCALARS'):
        data.append(line)
    f.close()
    data.sort()
    return data


  def _compare(self, file1, file2):
    def _convertLine(line):
      parts = line.split()
      return (int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]))
    data1 = self._readVTK(file1)
    data2 = self._readVTK(file2)
    if not data1 == data2:
      # Do full check
      ok = True
      for line1, line2 in zip(data1, data2):
        v1,x1,y1,z1 = _convertLine(line1)
        v2,x2,y2,z2 = _convertLine(line2)
        if not v1 == v2:
          ok = False
          print('ERROR: Nonmatching vertex sets')
        if abs(x1 - x2) > 1.0e-10:
          ok = False
          print('ERROR: Nonmatching x displacement, vertex %d, %g != %g' % \
                   (v1, x1, x2))
        if abs(y1 - y2) > 1.0e-10:
          ok = False
          print('ERROR: Nonmatching y displacement, vertex %d, %g != %g' % \
                   (v1, y1, y2))
        if abs(z1 - z2) > 1.0e-10:
          ok = False
          print('ERROR: Nonmatching z displacement, vertex %d, %g != %g' % \
                   (v1, z1, z2))
      if not ok:
        sys.exit("File '%s' DOES NOT MATCH file '%s'." % (file1, file2))
    print "File '%s' MATCHES '%s'." % (file1, file2)
    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = Verifier()
  app.run()

# End of file
