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
