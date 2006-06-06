#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2004  All Rights Reserved
#
#  Copyright 2004 Rensselaer Polytechnic Institute.
#  All worldwide rights reserved.  A license to use, copy, modify and
#  distribute this software for non-commercial research purposes only
#  is hereby granted, provided that this copyright notice and
#  accompanying disclaimer is not modified or removed from the software.
#
#  DISCLAIMER:  The software is distributed "AS IS" without any express
#  or implied warranty, including but not limited to, any implied
#  warranties of merchantability or fitness for a particular purpose
#  or any warranty of non-infringement of any current or pending patent
#  rights.  The authors of the software make no representations about
#  the suitability of this software for any particular purpose.  The
#  entire risk as to the quality and performance of the software is with
#  the user.  Should the software prove defective, the user assumes the
#  cost of all necessary servicing, repair or correction.  In
#  particular, neither Rensselaer Polytechnic Institute, nor the authors
#  of the software are liable for any indirect, special, consequential,
#  or incidental damages related to the software, to the maximum extent
#  the law permits.
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
        self.inventory.scanner.inventory.fileRoot, mesh = pylith3d.processMesh(self.inventory.scanner.inventory.fileRoot)
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
