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
# First feeble attempt at defining a Solver state.

from pyre.components.Component import Component

class SolverQSState(Component):
  """Python manager for state of implicit quasi-static time-stepping solver.
  The current setup assumes that all force and displacement vectors are
  PETSc Vecs, and the stiffness matrix is a PETSc Mat."""

  def setup(self, mesh, bc, solution):
    """Function to initialize components of Solver State.
    Present version uses the petscutils package, and it is assumed that
    PETSc is already initialized."""
    import pylith3d
    import petscutil.petscutil as petscbindings

    # At present, we are using Matt's function defined in scanner.cc.
    # It is unclear if we actually need separate vectors for rhs and sol,
    # but leave them in for now.
    stiff, rhs, sol = pylith3d.createPETScMat(mesh)

    # Set up solver state flags based on BC and solution method.
    if bc.numTractionBc != 0:
      self.solverStateFlags['tractionFlag'] = 1

    if solution.gravityX.value != 0.0 or \
         solution.gravityY.value != 0.0 or \
         solution.gravityZ.value != 0.0:
      self.solverStateFlags['gravityFlag'] = 1

    if bc.numConcForces != 0 or bc.numDiffForceEntries != 0:
      self.solverStateFlags['concForceFlag'] = 1

    if self.solverStateFlags['tractionFlag'] != 0 or \
         self.solverStateFlags['gravityFlag'] != 0 or \
         self.solverStateFlags['concForceFlag'] != 0:
      self.solverStateFlags['externFlag'] = 1
    
    if bc.numWinkForces != 0:
      self.solverStateFlags['winkFlag'] = 1

    if bc.numWinkxForces != 0:
      self.solverStateFlags['winkxFlag'] = 1

    if solution.usePreviousDisplacementFlag != 0:
      self.solverStateFlags['dprevFlag'] = 1

    # Allocate remaining PETSc vectors, depending on flag settings.
    self.solverState['petscStiff'] = stiff
    self.solverState['petscRhs'] = rhs
    self.solverState['petscSol'] = sol
    self.solverState['petscBintern'] = petscbindings.VecDuplicate(rhs)
    self.solverState['petscBresid'] = petscbindings.VecDuplicate(rhs)
    self.solverState['petscDispVec'] = petscbindings.VecDuplicate(rhs)

    if self.solverStateFlags['externFlag'] != 0:
      self.solverState['petscBextern'] = petscbindings.VecDuplicate(rhs)

    if self.solverStateFlags['tractionFlag'] != 0:
      self.solverState['petscBtraction'] = petscbindings.VecDuplicate(rhs)

    if self.solverStateFlags['gravityFlag'] != 0:
      self.solverState['petscBgravity'] = petscbindings.VecDuplicate(rhs)

    if self.solverStateFlags['concForceFlag'] != 0:
      self.solverState['petscBconcForce'] = petscbindings.VecDuplicate(rhs)

    if self.solverStateFlags['winkFlag'] != 0:
      self.solverState['petscBwink'] = petscbindings.VecDuplicate(rhs)

    if self.solverStateFlags['winkxFlag'] != 0:
      self.solverState['petscBwinkx'] = petscbindings.VecDuplicate(rhs)

    if self.solverStateFlags['dprevFlag'] != 0:
      self.solverState['petscDprev'] = petscbindings.VecDuplicate(rhs)

    return

  def assemble(self):
    """Assemble matrices and vectors."""

    import petscutil.petscutil as petscbindings

    import PETSc.Mat
    from PETSc.MatAssemblyType import MAT_FINAL_ASSEMBLY

    petscbindings.MatAssemblyBegin(self.solverState['petscStiff'],MAT_FINAL_ASSEMBLY)
    petscbindings.VecAssemblyBegin(self.solverState['petscRhs'])
    petscbindings.VecAssemblyBegin(self.solverState['petscSol'])
    if self.solverStateFlags['externFlag'] != 0:
      petscbindings.VecAssemblyBegin(self.solverState['petscBextern'])
    if self.solverStateFlags['tractionFlag'] != 0:
      petscbindings.VecAssemblyBegin(self.solverState['petscBtraction'])
    if self.solverStateFlags['gravityFlag'] != 0:
      petscbindings.VecAssemblyBegin(self.solverState['petscBgravity'])
    if self.solverStateFlags['concForceFlag'] != 0:
      petscbindings.VecAssemblyBegin(self.solverState['petscBconcForce'])
    petscbindings.VecAssemblyBegin(self.solverState['petscBintern'])
    petscbindings.VecAssemblyBegin(self.solverState['petscBresid'])
    if self.solverStateFlags['winkFlag'] != 0:
      petscbindings.VecAssemblyBegin(self.solverState['petscBwink'])
    if self.solverStateFlags['winkxFlag'] != 0:
      petscbindings.VecAssemblyBegin(self.solverState['petscBwinkx'])
    petscbindings.VecAssemblyBegin(self.solverState['petscDispVec'])
    if self.solverStateFlags['dprevFlag'] != 0:
      petscbindings.VecAssemblyBegin(self.solverState['petscDprev'])
    petscbindings.MatAssemblyEnd(self.solverState['petscStiff'],MAT_FINAL_ASSEMBLY)
    petscbindings.VecAssemblyEnd(self.solverState['petscRhs'])
    petscbindings.VecAssemblyEnd(self.solverState['petscSol'])
    if self.solverStateFlags['externFlag'] != 0:
      petscbindings.VecAssemblyEnd(self.solverState['petscBextern'])
    if self.solverStateFlags['tractionFlag'] != 0:
      petscbindings.VecAssemblyEnd(self.solverState['petscBtraction'])
    if self.solverStateFlags['gravityFlag'] != 0:
      petscbindings.VecAssemblyEnd(self.solverState['petscBgravity'])
    if self.solverStateFlags['concForceFlag'] != 0:
      petscbindings.VecAssemblyEnd(self.solverState['petscBconcForce'])
    petscbindings.VecAssemblyEnd(self.solverState['petscBintern'])
    petscbindings.VecAssemblyEnd(self.solverState['petscBresid'])
    if self.solverStateFlags['winkFlag'] != 0:
      petscbindings.VecAssemblyEnd(self.solverState['petscBwink'])
    if self.solverStateFlags['winkxFlag'] != 0:
      petscbindings.VecAssemblyEnd(self.solverState['petscBwinkx'])
    petscbindings.VecAssemblyEnd(self.solverState['petscDispVec'])
    if self.solverStateFlags['dprevFlag'] != 0:
      petscbindings.VecAssemblyEnd(self.solverState['petscDprev'])

    return

  def tearDown(self):
    """Destroy solver state vecs and mats."""

    import pylith3d
    import petscutil.petscutil as petscbindings

    pylith3d.destroyPETScMat(self.solverState['petscStiff'],
                               self.solverState['petscRhs'], 
                               self.solverState['petscSol'])
    if self.solverStateFlags['externFlag'] != 0:
      petscbindings.VecDestroy(self.solverState['petscBextern'])
    if self.solverStateFlags['tractionFlag'] != 0:
      petscbindings.VecDestroy(self.solverState['petscBtraction'])
    if self.solverStateFlags['gravityFlag'] != 0:
      petscbindings.VecDestroy(self.solverState['petscBgravity'])
    if self.solverStateFlags['concForceFlag'] != 0:
      petscbindings.VecDestroy(self.solverState['petscBconcForce'])
    petscbindings.VecDestroy(self.solverState['petscBintern'])
    petscbindings.VecDestroy(self.solverState['petscBresid'])
    if self.solverStateFlags['winkFlag'] != 0:
      petscbindings.VecDestroy(self.solverState['petscBwink'])
    if self.solverStateFlags['winkxFlag'] != 0:
      petscbindings.VecDestroy(self.solverState['petscBwinkx'])
    petscbindings.VecDestroy(self.solverState['petscDispVec'])
    if self.solverStateFlags['dprevFlag'] != 0:
      petscbindings.VecDestroy(self.solverState['petscDprev'])

    return
    
  def __init__(self):
    """Constructor."""
    # Dictionary containing all Solver State objects.
    self.solverState = {'petscStiff': None,
                        'petscRhs': None,
                        'petscSol': None,
                        'petscBextern': None,
                        'petscBtraction': None,
                        'petscBgravity': None,
                        'petscBconcForce': None,
                        'petscBintern': None,
                        'petscBresid': None,
                        'petscBwink': None,
                        'petscBwinkx': None,
                        'petscDispVec': None,
                        'petscDprev': None}

    self.solverStateFlags = {'externFlag': 0,
                             'tractionFlag': 0,
                             'gravityFlag': 0,
                             'concForceFlag': 0,
                             'winkFlag': 0,
                             'winkxFlag': 0,
                             'dprevFlag': 0}
    return
        

# version
__id__ = "$SolverQSState.py$"

# End of file 
