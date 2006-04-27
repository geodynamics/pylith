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

# Local base class for Pyre simulation solver derived from pyre SimpleSolver class.
# This version is based almost entirely on Brad Aagaard's base Solver class for
# EqSim.

from pyre.simulations.SimpleSolver import SimpleSolver as SolverBase

class Solver(SolverBase):
  """Base class for pylith simulation solver."""

  # INVENTORY

  class Inventory(SolverBase.Inventory):
    """Python object for managing Solver facilities and properties."""

    import pyre.inventory

    from Modeler import Modeler
    modeler = pyre.inventory.facility("modeler", factory=Modeler)
    modeler.meta['tip'] = "3D solid geometry."

    from Mesher import Mesher
    mesher = pyre.inventory.facility("mesher", factory=Mesher)
    mesher.meta['tip'] = "Finite element mesh generator or importer."

    # stuff dependent on petscutil may change
    from petscutil.PetscManager import PetscManager
    petsc = pyre.inventory.facility("petsc", factory=PetscManager)
    petsc.meta['tip'] = "Manager of PETSc routines/flags."

    from petscutil.PetscLogger import PetscLogger
    logger = pyre.inventory.facility("logger", factory=PetscLogger)
    logger.meta['tip'] = "Manager of PETSc logging."

    # PUBLIC METHODS

    def launch(self, app):
      """Launch solver."""
      SolverBase.launch(self, app)

      # Get communicator
      communicator = app.layout.communicator
      self._petsc.setWorldComm(communicator)
      self._petsc.initialize(self.name)
      return

    def endTimestep(self, t):
      """Complete time step."""
      SolverBase.endTimestep(self, t)
      return

    def endSimulation(self, steps, t):
      """Clean up at end of simulation."""
      SolverBase.endSimulation(self, steps, t)
      self._petsc.finalize()
      return

    def __init__(self, name="solver"):
      """Constructor."""
      SolverBase.__init__(self, name)
      self.state = None
      self.mesh = None
      self.modeler = None
      self.mesher = None
      return


    # PRIVATE METHODS

    def _configure(self):
      """Setup members using inventory."""
      self._petsc = self.inventory.petsc
      self._logger = self.inventory.logger
      self.modeler = self.inventory.modeler
      self.mesher = self.inventory.mesher
      return

# version
__id__ = "$Id$"

# End of file 
