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


from pyre.applications.Application import Application as BaseApplication

import journal

class Application(BaseApplication):


    def run(self):
        self._progress.log("running application %s" % self.name)
#        from time import clock as now
#        start = now()
        self.initialize()
        self.launch()
#        finish = now()
#        usertime = finish - start
#        print "Total user time:  %g" % usertime
        return

    def initialize(self):
        self._progress.log("initializing application %s" % self.name)
        self.layout()
        lm3dsetup = self.inventory.setup
        lm3dsetup.initialize(self.inventory.scanner)
        lm3dsetup.run()
        lm3drun = self.inventory.solverold
        lm3drun.initialize(self.inventory.scanner, self.inventory.setup)
        # This needs to be changed to be part of Solver.
        # lm3drun.run()
        return

    def launch(self):
        self._progress.log("launching application %s" % self.name)
        controller = self.inventory.controller
        controller.launch(self)
        totalTime=self.inventory.total_time
        controller.march(totalTime)
        return

    def layout(self):
        self._progress.log("%s layout" % self.name)
        layout = self.inventory.layout
        layout.layout(self)
        self.solver = self.inventory.solver
        return

    def reportConfiguration(self):
        self._info.line("configuration:")
        self._info.line("facilities:")
        self._info.line("  layout: %s" % self.inventory.layout)
        self._info.line("  controller: %s" % self.inventory.controller)
        self._info.line("  solver: %s" % self.inventory.solver)
        self._info.line("  total_time: %s" % self.inventory.total_time)
        self._info.log()
        return


    def __init__(self, name="lithomop3d"):
        BaseApplication.__init__(self, name)
        self.solver = None

        self._progress = journal.debug("%s.progress" % name)
        self._info = journal.info("%s.info" % name)
        return


    class Inventory(BaseApplication.Inventory):

        import pyre.facilities
        import pyre.properties
        from LayoutManager import LayoutManager
        from Controller import Controller
        from Solver import Solver
        from pyre.units.time import second
        from Lithomop3d_scan import Lithomop3d_scan
        from Lithomop3d_setup import Lithomop3d_setup
        from Lithomop3d_run import Lithomop3d_run

        inventory = [
            pyre.facilities.facility("scanner", default=Lithomop3d_scan()),
            pyre.facilities.facility("setup", default=Lithomop3d_setup()),
            pyre.facilities.facility("solverold", default=Lithomop3d_run()),
            pyre.facilities.facility("layout", default=LayoutManager()),
            pyre.facilities.controller(default=Controller()),
            pyre.facilities.solver("solver",default=Solver()),
            pyre.properties.dimensional("total_time",default=0.0*second),
            ]


# version
# $Id: Application.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $

# End of file 
