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

## @file verification/verify

## @brief Python application to check VTK files from 1 processor and 2
## processor runs of PyLith.

import os, re, sys

from pyre.applications.Script import Script as Application

class Verifier(Application):
  """
  Python application to check VTK files from 1 processor and 2
  processors runs of PyLith.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing Verifier facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Verifier facilities and properties.
    ##
    ## \b Properties
    ## @li \b cfg_file Filename of .cfg file with simulation parameters.
    ## @li \b pylith Absolute path to 'pylith' executable.
    ## @li \b petsc_options PETSc options.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    cfgFile = pyre.inventory.str("cfg_file", default="pylithapp.cfg")
    cfgFile.meta['tip'] = "Filename of .cfg file with simulation parameters."

    pylith = pyre.inventory.str("pylith", default="`which pylith`")
    pylith.meta['tip'] = "Absolute path to 'pylith' executable."

    petscOptions = pyre.inventory.list("petsc_options",
                                        default=['--petsc.ksp_type=preonly',
                                                 '--petsc.pc_type=lu',
                                                 '--petsc.mat_type=aijmumps'])
    petscOptions.meta['tip'] = "PETSc options."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="verifyapp"):
    Application.__init__(self, name)
    self.exp = re.compile(r'SCALARS \w+_verify_t. double 4')
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
