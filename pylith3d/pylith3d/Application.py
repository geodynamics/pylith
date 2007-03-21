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


from cig.cs.petsc import PetscApplication


class Application(PetscApplication):


    name = "pylith3d"


    # Tell the framework where to find PETSc functions.
    import pylith3d as petsc


    # Use PETSc-style command line parsing.
    from cig.cs.petsc import PetscCommandlineParser as CommandlineParser


    def main(self, *args, **kwds):
    
#        from time import clock as now
#        start = now()
        pl3dsetup = self.inventory.setup
        import pylith3d

        scanner = self.inventory.scanner

        mesh = pylith3d.processMesh(scanner.macroString(scanner.Inventory.outputFileRoot),
                                    scanner.macroString(scanner.Inventory.inputFileRoot),
                                    scanner.inventory.interpolateMesh,
                                    scanner.inventory.partitioner)

        scanner.initialize()
        
        pl3dsetup.initialize(scanner)
        pl3dsetup.read()
        pl3dsetup.numberequations()
        pl3dsetup.sortmesh()
        pl3dsetup.sparsesetup(mesh)
        pl3dsetup.allocateremaining()
        pl3dsetup.meshwrite()
        pl3drun = self.inventory.solver
        pl3drun.fileRoot = scanner.inventory.fileRoot
        pl3drun.pointerToIelindx = pl3dsetup.pointerToIelindx
        pl3drun.mesh = mesh
        pl3drun.initialize(scanner, self.inventory.setup)
        pl3drun.run()
#        finish = now()
#        usertime = finish - start
#        print "Total user time:  %g" % usertime
        return


    class Inventory(PetscApplication.Inventory):

        import pyre.inventory
        from cig.cs.petsc import PetscProperty
        from Pylith3d_scan import Pylith3d_scan
        from Pylith3d_setup import Pylith3d_setup
        from Pylith3d_run import Pylith3d_run

        scanner = pyre.inventory.facility("scanner", factory=Pylith3d_scan)
        setup = pyre.inventory.facility("setup", factory=Pylith3d_setup)
        solver = pyre.inventory.facility("solver", factory=Pylith3d_run)

        # declare PETSc options that are of interest to PyLith
        ksp_monitor        = PetscProperty()
        ksp_view           = PetscProperty()
        ksp_rtol           = PetscProperty()
        log_summary        = PetscProperty()
        pc_type            = PetscProperty()
        sub_pc_type        = PetscProperty()
        start_in_debugger  = PetscProperty()
        debugger_pause     = PetscProperty()


# version
# $Id: Application.py,v 1.5 2005/04/15 00:18:21 willic3 Exp $

# End of file 
