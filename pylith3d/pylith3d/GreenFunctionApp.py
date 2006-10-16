#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
from pyre.applications.Script import Script as BaseScript

class GreenFunctionApp(BaseScript):
    def readSamplePoints(self, filename):
        '''Read in the sampling locations
        - One point per line, three values per line (x,y,z)
        - Returns a Numeric array'''
        import Numeric
        f = file(filename)
        points = []
        for line in f.readlines():
            points.append([float(v) for v in line.strip().split(' ')])
        f.close()
        return Numeric.array(points)

    def outputSampleValues(self, filename, values):
        '''sample# sample values impluse# impulse type'''
        # Computing normal to the fault:
        #   Split nodes define the fault
        #   Get all fault faces for a node
        #   Area weighted average of normals
        f = file(filename, 'w')
        for v, values in enumerate(values):
            write(f, '%d %g %g %g 1 0' % (v, values[0], values[1], values[2]))
        f.close()
        return

    def main(self, *args, **kwds):
        import pylith3d

        pl3dscanner = self.inventory.scanner
        pl3dsetup   = self.inventory.setup
        pl3drun     = self.inventory.solver
        points      = readSamplePoints()
        pylith3d.PetscInitialize(pl3dscanner.inventory.fileRoot+'.sample')
        self.inventory.scanner.inventory.fileRoot, mesh = pylith3d.processMesh(pl3dscanner.inventory.fileRoot, pl3dscanner.inventory.interpolateMesh)
        try:
            pl3dsetup.initialize(pl3dscanner)
        except self.inventory.scanner.CanNotOpenInputOutputFilesError, error:
            import sys
            print >> sys.stderr
            error.report(sys.stderr)
            print >> sys.stderr
            print >> sys.stderr, "%s: %s" % (error.__class__.__name__, error)
            sys.exit(1)
        pl3dsetup.read()
        pl3dsetup.numberequations()
        pl3dsetup.sortmesh()
        pl3dsetup.sparsesetup(mesh)
        pl3dsetup.allocateremaining()
        pl3dsetup.meshwrite()
        pl3drun.fileRoot = self.inventory.scanner.inventory.fileRoot
        pl3drun.pointerToIelindx = pl3dsetup.pointerToIelindx
        pl3drun.mesh = mesh
        pl3drun.initialize(self.inventory.scanner, self.inventory.setup)
        pl3drun.solveElastic()
        values = pl3drun.interpolatePoints(points)
        self.outputSampleValues(pl3dscanner.inventory.fileRoot+'.output', values):
        return


    def __init__(self, name = "pylith3d"):
        BaseScript.__init__(self, name)
        return

    class Inventory(BaseScript.Inventory):
        import pyre.inventory
        from Pylith3d_scan  import Pylith3d_scan
        from Pylith3d_setup import Pylith3d_setup
        from Pylith3d_run   import Pylith3d_run

        scanner = pyre.inventory.facility("scanner", factory = Pylith3d_scan)
        setup   = pyre.inventory.facility("setup",   factory = Pylith3d_setup)
        solver  = pyre.inventory.facility("solver",  factory = Pylith3d_run)

# version
# $Id: Application.py,v 1.5 2005/04/15 00:18:21 willic3 Exp $

# End of file 
