#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#             Copyright (C) 2006 Rensselaer Polytechnic Institute
#
# 
# 	Permission is hereby granted, free of charge, to any person
# 	obtaining a copy of this software and associated documentation
# 	files (the "Software"), to deal in the Software without
# 	restriction, including without limitation the rights to use,
# 	copy, modify, merge, publish, distribute, sublicense, and/or
# 	sell copies of the Software, and to permit persons to whom the
# 	Software is furnished to do so, subject to the following
# 	conditions:
# 
# 	The above copyright notice and this permission notice shall be
# 	included in all copies or substantial portions of the Software.
# 
# 	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# 	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# 	OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# 	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# 	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# 	WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# 	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# 	OTHER DEALINGS IN THE SOFTWARE.
#         
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# First feeble attempt at defining a Solver state.

from pyre.components.Component import Component

class SolverQSState(Component):
  """Python manager for state of implicit quasi-static time-stepping solver."""

  def __init__(self):
    """Constructor."""
    # Create dictionary for force info.
    self.force = {'sizeBextern': 0,
                  'ptrBextern': None,
                  'sizeBtraction': 0,
                  'ptrBtraction': None,
                  'sizeBgravity': 0,
                  'ptrBgravity': None,
                  'sizeBconcForce': 0,
                  'ptrBconcForce': None,
                  'sizeBintern': 0,
                  'ptrBintern': None,
                  'sizeBresid': 0,
                  'ptrBresid': None,
                  'sizeBwink': 0,
                  'ptrBwink': None,
                  'sizeBwinkx': 0,
                  'ptrBwinkx': None}

    # Create dictionary for displacement vector info.
    self.dispvec = {'sizeDispVec': 0,
                    'ptrDispVec': None,
                    'sizeDprev': 0,
                    'ptrDprev': None}

    # Create dictionary for stiffness, rhs, and sol.
    # I need to figure out how this should really work.  Matt is defining a single
    # PETSc 'mesh' container that includes A (global stiffness), rhs (PETSc Vec), and
    # sol (PETSc Vec).  It is quite likely that all of the vectors above will be
    # gathered into a single 'solverState' container.
    self.stiff = {'petscStiff': None,
                  'petscRhs': None,
                  'petscSol': None}

    return
        

# version
__id__ = "$Id$"

# End of file 
