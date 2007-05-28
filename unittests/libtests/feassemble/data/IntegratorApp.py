#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/IntegratorApp.py

## @brief Python application for generating C++ data files for testing
## C++ integrator objects.

from pyre.applications.Script import Script

import numpy

# IntegratorApp class
class IntegratorApp(Script):
  """
  Python application for generating C++ data files for testing C++
  integrator objects.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Script.Inventory):
    """Python object for managing IntegratorApp facilities and properties."""

    ## @class Inventory
    ## Python object for managing IntegratorApp facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b data Data manager.

    import pyre.inventory

    from pylith.utils.CppData import CppData
    data = pyre.inventory.facility("data", factory=CppData)
    data.meta['tip'] = "Data manager."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="integratorapp"):
    """
    Constructor.
    """
    Script.__init__(self, name)

    # Mesh information
    self.meshFilename = None

    # Quadrature information
    self.spaceDim = None
    self.cellDim = None
    self.numBasis = None
    self.numQuadPts = None
    self.quadPts = None
    self.quadWts = None
    self.basis = None
    self.basisDeriv = None

    # Material information
    self.matType = None
    self.matDBFilename = None
    self.matId = None
    self.matLabel = None

    # Input fields
    self.dt = None
    self.fieldTpdt = None
    self.fieldT = None
    self.fieldTmdt = None

    # Calculated values
    self.valsResidual = None
    self.valsJacobian = None
    return


  def main(self):
    """
    Run the application.
    """
    self._initialize()
    self._calculateResidual()
    self._calculateJacobian()
    self._initData()
    self.data.write(self.name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members using inventory.
    """
    Script._configure(self)
    self.data = self.inventory.data
    return


  def _initialize(self):
    """
    Get quadrature information.
    """
    q = self.quadrature
    q.calculateBasis()
    self.spaceDim = q.spaceDim
    self.numBasis = q.numBasis
    self.cellDim = q.cellDim
    self.numQuadPts = q.numQuadPts
    self.basis = q.basis
    self.basisDeriv = q.basisDeriv
    self.quadPts = q.quadPtsRef
    self.quadWts = q.quadWts

    return
  

  def _calculateResidual(self):
    """
    Calculate contribution to residual of operator for integrator.
    """
    raise NotImplementedError("Not implemented in abstract base class.")


  def _calculateJacobian(self):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.
    """
    self.valsAction = numpy.dot(self.valsMatrix, self.fieldIn)[:]
    return
    

  def _initData(self):
    # Mesh information
    self.data.addScalar(vtype="char*", name="_meshFilename",
                        value=self.meshFilename,
                        format="%s")

    # Quadrature information
    self.data.addScalar(vtype="int", name="_spaceDim", value=self.spaceDim,
                        format="%d")
    self.data.addScalar(vtype="int", name="_cellDim", value=self.cellDim,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numBasis", value=self.numBasis,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numQuadPts", value=self.numQuadPts,
                        format="%d")
    self.data.addArray(vtype="double", name="_quadPts", values=self.quadPts,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_quadWts", values=self.quadWts,
                       format="%16.8e", ncols=self.numQuadPts)
    
    self.data.addArray(vtype="double", name="_basis", values=self.basis,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_basisDeriv",
                       values=self.basisDeriv,
                       format="%16.8e", ncols=self.cellDim)

    # Material information
    self.data.addScalar(vtype="char*", name="_matType", value=self.matType,
                        format="%s")
    self.data.addScalar(vtype="char*", name="_matDBFilename",
                        value=self.matDBFilename, format="%s")
    self.data.addScalar(vtype="int", name="_matId",
                        value=self.matId, format="%d")
    self.data.addScalar(vtype="char*", name="_matLabel",
                        value=self.matLabel, format="%s")

    # Input files
    self.data.addScalar(vtype="double", name="_dt", values=self.dt,
                        format="%16.8e")
    self.data.addArray(vtype="double", name="_fieldTpdt",
                       values=self.fieldTpdt,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="double", name="_fieldT",
                       values=self.fieldTpdt,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="double", name="_fieldTmdt",
                       values=self.fieldTpdt,
                       format="%16.8e", ncols=self.spaceDim)

    # Calculated values
    self.data.addArray(vtype="double", name="_valsResidual",
                       values=self.valsResidual,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="double", name="_valsJacobian",
                       values=self.valsJacobian,
                       format="%16.8e", ncols=self.spaceDim)
      
    return

  
# End of file 
