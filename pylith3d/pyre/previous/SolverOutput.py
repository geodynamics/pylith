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
