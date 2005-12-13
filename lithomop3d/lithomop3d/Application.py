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
        lm3dsetup = self.inventory.setup
        import lithomop3d
        lithomop3d.PetscInitialize()
        self.inventory.scanner.inventory.fileRoot, mesh = lithomop3d.processMesh(self.inventory.scanner.inventory.fileRoot)
        lm3dsetup.initialize(self.inventory.scanner)
        lm3dsetup.read()
        lm3dsetup.numberequations()
        lm3dsetup.sortmesh()
        lm3dsetup.sparsesetup(mesh)
        lm3dsetup.allocateremaining()
        lm3dsetup.meshwrite()
        lm3drun = self.inventory.solver
        lm3drun.initialize(self.inventory.scanner, self.inventory.setup)
        lm3drun.run()
#        finish = now()
#        usertime = finish - start
#        print "Total user time:  %g" % usertime
        return


    def __init__(self, name="lithomop3d"):
        BaseScript.__init__(self, name)
        return


    class Inventory(BaseScript.Inventory):

        import pyre.inventory
        from Lithomop3d_scan import Lithomop3d_scan
        from Lithomop3d_setup import Lithomop3d_setup
        from Lithomop3d_run import Lithomop3d_run

        scanner = pyre.inventory.facility("scanner", factory=Lithomop3d_scan)
        setup = pyre.inventory.facility("setup", factory=Lithomop3d_setup)
        solver = pyre.inventory.facility("solver", factory=Lithomop3d_run)
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
