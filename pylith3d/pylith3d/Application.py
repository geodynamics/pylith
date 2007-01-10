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


from mpi.Application import Application as BaseApplication


class Application(BaseApplication):


    def main(self, *args, **kwds):
        import sys
        
        # if we are embedding, insert the extension module in the
        # 'pylith3d' package
        try:
            import builtin_pylith3d
            sys.modules['pylith3d.pylith3d'] = builtin_pylith3d
        except ImportError:
            pass
    
#        from time import clock as now
#        start = now()
        pl3dsetup = self.inventory.setup
        import pylith3d
        pylith3d.PetscInitialize(self.petscArgs)
        self.inventory.scanner.inventory.fileRoot, mesh = pylith3d.processMesh(self.inventory.scanner.inventory.fileRoot, self.inventory.scanner.inventory.interpolateMesh, self.inventory.scanner.inventory.partitioner)
        try:
            pl3dsetup.initialize(self.inventory.scanner)
        except self.inventory.scanner.CanNotOpenInputOutputFilesError, error:
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
        BaseApplication.__init__(self, name)
        return


    def _configure(self):
        self.petscArgs = []
        
        ksp_monitor = self.inventory.ksp_monitor
        if ksp_monitor:
            self.petscArgs.append("-ksp_monitor")
            if ksp_monitor != "true":
                self.petscArgs.append(ksp_monitor)
        
        if self.inventory.ksp_view:
            self.petscArgs.append("-ksp_view")
        
        ksp_rtol = self.inventory.ksp_rtol
        if ksp_rtol:
            self.petscArgs.extend(["-ksp_rtol", ksp_rtol])
        
        if self.inventory.log_summary:
            self.petscArgs.append("-log_summary")

        pc_type = self.inventory.pc_type
        sub_pc_type = self.inventory.sub_pc_type
        self.petscArgs.extend(["-pc_type", pc_type, "-sub_pc_type", sub_pc_type])

        if self.inventory.start_in_debugger:
            self.petscArgs.append("-start_in_debugger")

        self.petscArgs.extend(self.inventory.petsc.getArgs())

        return


    def createCommandlineParser(self):
        from PetscUtil import PetscCommandlineParser
        return PetscCommandlineParser()


    class Inventory(BaseApplication.Inventory):

        import pyre.inventory
        from Pylith3d_scan import Pylith3d_scan
        from Pylith3d_setup import Pylith3d_setup
        from Pylith3d_run import Pylith3d_run
        from PetscUtil import PetscFacility

        scanner = pyre.inventory.facility("scanner", factory=Pylith3d_scan)
        setup = pyre.inventory.facility("setup", factory=Pylith3d_setup)
        solver = pyre.inventory.facility("solver", factory=Pylith3d_run)

        # declare PETSc options that are of interest to PyLith
        ksp_monitor = pyre.inventory.str("ksp_monitor")
        ksp_view = pyre.inventory.bool("ksp_view")
        ksp_rtol = pyre.inventory.str("ksp_rtol")
        log_summary = pyre.inventory.bool("log_summary")
        pc_type = pyre.inventory.str("pc_type")
        sub_pc_type = pyre.inventory.str("sub_pc_type")
        start_in_debugger = pyre.inventory.str("start_in_debugger")

        # a dummy facility for passing arbitrary options to PETSc
        petsc = PetscFacility("petsc")


# version
# $Id: Application.py,v 1.5 2005/04/15 00:18:21 willic3 Exp $

# End of file 
