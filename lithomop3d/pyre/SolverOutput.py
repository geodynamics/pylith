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

from pyre.components.Component import Component

import journal

# ----------------------------------------------------------------------
class SolverOutput(Component):

  def outputNow(self, t):

    if self._prevT == None:
      outputFlag = True
    else:
      timeStep = self.inventory.outputTimeStep
      outputFlag = t-self._prevT > timeStep
    return outputFlag

  def writeTimeStep(self, t, XXXXXXX):
    if self.outputnow(t):
      self._writeTimeStep(t, XXXXXXX)
      self._prevT = t
    return

  def _writeTimeStep(self, t, XXXXXXX):
    return
      
  def __init__(self):

    Component.__init__(self, "SolverOutput", "SolverOutput")

    self._progress = journal.debug("SolverOutput.progress")
    self._prevT = None
    return

  class Inventory(Component.Inventory):
      
    import pyre.properties as props
    from pyre.units.time import year
    
    inventory = [
      props.dimensional("outputTimeStep", default=10.0*year)
      ]

# version
__id__ = "$Id: SolverOutput.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $"

# End of file 
