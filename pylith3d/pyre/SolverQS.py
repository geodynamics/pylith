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

from Solver import Solver as SolverBase


class SolverQS(SolverBase):
  """Quasi-static solver."""

  # INVENTORY
  
  class Inventory(SolverBase.Inventory):
    """Python object for managing Solver facilities and properties."""

    import pyre.inventory

    # No inventory yet.

  def launch(self, app):
    """Launch solver."""
    SolverBase.launch(self, app)

    # Setup PETSc logging
    self._logger.registerStages(['launch', 'newStep',
                                 'advance', 'endTimestep'])

    self._logger.registerEvents(['iterate'])
    # Need to add more events, and make sure this method will work correctly.
    # Also, I will probably eventually switch over to PETSc's nonlinear solver,
    # which will replace iterate.
    
    self._logger.stagePush('launch')

    # Setup the model infrastructure
    self._initializeModel()

    # Setup the solver state
    self._initializeState()

    self._logger.stagePop()
    return

  # Need to start adding in functions to replace my current fortran.
  def newStep(self, t, step):
    """Presently set up to zero state (and state increment) variables
    at the beginning of a time step using PETSc VecScale routine."""
    self._logger.stagePush('newStep')

    import petscutil.petscutil as petscbindings
    petscbindings.VecScale(

    
    return


  def __init__(self, name="solverqs"):
    """Constructor."""
    SolverBase.__init__(self, name)

    from SolverQSState import SolverQSState
    self.state = SolverQSState()

    import journal
    self._info = journal.info(name)
        
    return


# version
__id__ = "$Id: SolverQS.py,v 1.1.1.1 2005/03/08 16:13:46 aivazis Exp $"

#
# End of file
