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


class Application(BaseScript):


    def main(self, *args, **kwds):
#        from time import clock as now
#        start = now()
        pl3dsetup = self.inventory.setup
        import pylith3d
        pylith3d.PetscInitialize()
        self.inventory.scanner.inventory.fileRoot, mesh = pylith3d.processMesh(self.inventory.scanner.inventory.fileRoot, self.inventory.scanner.inventory.interpolateMesh, self.inventory.scanner.inventory.partitioner)
        try:
            pl3dsetup.initialize(self.inventory.scanner)
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
        pl3drun = self.inventory.solver
        pl3drun.fileRoot = self.inventory.scanner.inventory.fileRoot
        pl3drun.pointerToIelindx = pl3dsetup.pointerToIelindx
        pl3drun.mesh = mesh
        pl3drun.initialize(self.inventory.scanner, self.inventory.setup)
        pl3drun.run()
#        finish = now()
#        usertime = finish - start
#        print "Total user time:  %g" % usertime
        return


    def __init__(self, name="pylith3d"):
        BaseScript.__init__(self, name)
        return


    class Inventory(BaseScript.Inventory):

        import pyre.inventory
        from Pylith3d_scan import Pylith3d_scan
        from Pylith3d_setup import Pylith3d_setup
        from Pylith3d_run import Pylith3d_run

        scanner = pyre.inventory.facility("scanner", factory=Pylith3d_scan)
        setup = pyre.inventory.facility("setup", factory=Pylith3d_setup)
        solver = pyre.inventory.facility("solver", factory=Pylith3d_run)
        # ksp_monitor = pyre.inventory.str("ksp_monitor",default="1")
        # ksp_view = pyre.inventory.str("ksp_view",default="1")
        # log_summary = pyre.inventory.str("log_summary",default="1")
        # log_info = pyre.inventory.str("log_info",default="0")
        # pc_type = pyre.inventory.str("pc_type",default="1")
        # start_in_debugger = pyre.inventory.str("start_in_debugger",default="0")
        # petsc_solver = pyre.inventory.str("petsc_solver",default="1")


# version
# $Id: Application.py,v 1.5 2005/04/15 00:18:21 willic3 Exp $

# End of file 
