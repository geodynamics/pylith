#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams
#  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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


from pyre.applications.Application import Application as BaseApplication


class Application(BaseApplication):


    def run(self, *args, **kwds):
#        from time import clock as now
#        start = now()
        pl3dsetup = self.inventory.setup
        pl3dsetup.initialize(self.inventory.scanner)
        pl3dsetup.run()
        pl3drun = self.inventory.solver
        pl3drun.initialize(self.inventory.scanner, self.inventory.setup)
        pl3drun.run()
#        finish = now()
#        usertime = finish - start
#        print "Total user time:  %g" % usertime
        return


    def __init__(self, name="pylith3d"):
        BaseApplication.__init__(self, name)
        return


    class Inventory(BaseApplication.Inventory):

        import pyre.facilities
        from Pylith3d_scan import Pylith3d_scan
        from Pylith3d_setup import Pylith3d_setup
        from Pylith3d_run import Pylith3d_run

        inventory = [
            pyre.facilities.facility("scanner", default=Pylith3d_scan()),
            pyre.facilities.facility("setup", default=Pylith3d_setup()),
            pyre.facilities.facility("solver", default=Pylith3d_run())
            ]


# version
# $Id: Application-old.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $

# End of file 
