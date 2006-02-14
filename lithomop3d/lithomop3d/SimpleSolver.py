#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2005  All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

from Solver import Solver


class SimpleSolver(Solver):


    class Inventory(Solver.Inventory):

        import pyre.inventory
        from pyre.units.SI import second

        timestep = pyre.inventory.dimensional("timestep", default=1.0e-6 * second)


    def stableTimestep(self):
        dt = self.inventory.timestep
        Solver.stableTimestep(self, dt)
        return dt


    def __init__(self, name=None):
        if name is None:
            name = "simpleSolver"

        Solver.__init__(self, name)
        return


# version
__id__ = "$Id: SimpleSolver.py,v 1.1.1.1 2005/03/08 16:13:46 aivazis Exp $"

#
# End of file
