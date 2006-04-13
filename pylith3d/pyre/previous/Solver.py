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

from pyre.components.SimpleSolver import SimpleSolver

import journal

import pylith3d as bindings

# ----------------------------------------------------------------
class Solver(SimpleSolver):

  def launch(self, application):
    SimpleSolver.launch(self, application)

    # This part will look a lot like Pylith3d_run.py.
    # It will be much easier to bundle things into python objects and then
    # unbundle them in the bindings.
    bindings.solver_const(XXXXXXXXXX)
    bindings.solver_elastc(XXXXXXXXXX)
    
    from pyre.units.time import second
    t, step = 0.0*second, 0
    self.save(t, XXXXXX)
    return (t, step)

  def newStep(self, t, step):
    SimpleSolver.newStep(self, t, step)
    bindings.solver_starttimestep(XXXXXXXXXXXXXXX)
    return

  def applyBoundaryConditions(self):
    SimpleSolver.applyBoundaryConditions(self)
    bindings.solver_applybc(XXXXXXXXXXXXX)
    return

  def advance(self, dt):
    SimpleSolver.advance(self, dt)
    bindings.solver_advance(XXXXXXXXXXXXX)
    return

  def save(self, t):
    output = self.inventory.output
    output.writeTimeStep(t, XXXXXXXXXXX)
    return
    
  def endTimestep(self, t):
    SimpleSolver.endTimestep(self, t)
    return

  def __init__(self):
    SimpleSolver.__init__(self, "Solver", "Solver")
    return

  class Inventory(SimpleSolver.Inventory):

    import pyre.facilities
    from SolverOutputAsciiUcd import SolverOutputAsciiUcd
    
    inventory = [
      pyre.facilities.facility("output", default=SolverOutputAsciiUcd())
      ]
    
# version
__id__ = "$Id: Solver.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $"

# End of file 
