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
