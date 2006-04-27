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

from pyre.components.SimulationController import SimulationController

class Controller(SimulationController):

  def launch(self, application):
    SimulationController.launch(self, application)
    self.solver = application.solver
    self.clock, self.step = self.solver.launch(application)
    return

  def march(self, totalTime=0, steps=0):
    """explicit time loop"""

    # the main simulation time loop
    while 1:
      
      # notify solvers we are starting a new timestep
      self.startTimestep()
      
      # compute an acceptable timestep
      dt = self.stableTimestep()
      
      # synchronize boundary information
      self.applyBoundaryConditions()
      
      # advance
      self.advance(dt)
      
      # update smulation clock and step number
      self.clock = self.clock + dt
      self.step = self.step + 1
      
      # do io
      self.save()

      # notify solver we finished a timestep
      self.endTimestep()
      
      # are we done?
      if totalTime and self.clock >= totalTime:
        break
      if steps and self.step >= steps:
        break
      
    # end of time advance loop           
      
    # Notify solver we are done
    self.endSimulation()
      
    return

  def save(self):
    self.solver.save(self.clock)
    return

  def __init__(self, name="Controller", requirements=None):
    SimulationController.__init__(self, name, requirements)
    return

  class Inventory(SimulationController.Inventory):

    inventory = [
      ]

# version
__id__ = "$Id: Controller.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $"

# End of file 
