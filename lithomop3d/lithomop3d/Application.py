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


class Application(BaseApplication):


    def run(self, *args, **kwds):
#        from time import clock as now
#        start = now()
        lm3dsetup = self.inventory.setup
        lm3dsetup.initialize(self.inventory.scanner)
        lm3dsetup.run()
        lm3drun = self.inventory.solver
        lm3drun.initialize(self.inventory.scanner, self.inventory.setup)
        lm3drun.run()
#        finish = now()
#        usertime = finish - start
#        print "Total user time:  %g" % usertime
        return


    def __init__(self, name="lithomop3d"):
        BaseApplication.__init__(self, name)
        return


    class Inventory(BaseApplication.Inventory):

        import pyre.facilities
        from Lithomop3d_scan import Lithomop3d_scan
        from Lithomop3d_setup import Lithomop3d_setup
        from Lithomop3d_run import Lithomop3d_run

        inventory = [
            pyre.facilities.facility("scanner", default=Lithomop3d_scan()),
            pyre.facilities.facility("setup", default=Lithomop3d_setup()),
            pyre.facilities.facility("solver", default=Lithomop3d_run()),
            pyre.properties.str("ksp_monitor",default="1"),
            pyre.properties.str("ksp_view",default="1"),
            pyre.properties.str("log_summary",default="1"),
            pyre.properties.str("log_info",default="0"),
            pyre.properties.str("pc_type",default="jacobi"),
            pyre.properties.str("start_in_debugger",default="0"),
            pyre.properties.str("petsc_solver",default="0"),
            ]


# version
# $Id: Application.py,v 1.3 2005/03/10 01:10:37 knepley Exp $

# End of file 
