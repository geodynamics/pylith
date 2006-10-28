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

    self.numVertices = None
    self.spaceDim = None
    self.numCells = None
    self.cellDim = None
    self.numCorners = None
    self.numQuadPts = None
    self.fiberDim = None

    self.quadPts = None
    self.quadWts = None
    self.vertices = None
    self.cells = None

    self.basis = None
    self.basisDeriv = None
    self.fieldIn = None
    self.valsActions = None
    self.valsMatrix = None
    return


  def main(self):
    """
    Run the application.
    """
    self._initialize()
    self._calculateMatrix()
    self._calculateAction()
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
    self.numCorners = q.numCorners
    self.cellDim = q.cellDim
    self.numQuadPts = q.numQuadPts
    self.basis = q.basis
    self.basisDeriv = q.basisDeriv
    self.quadPts = q.quadPtsRef
    self.quadWts = q.quadWts

    return
  

  def _calculateMatrix(self):
    """
    Calculate matrix associated with integration.
    """
    raise NotImplementedError("Not implemented in abstract base class.")


  def _calculateAction(self):
    """
    Calculate integration action using matrix.
    """
    self.valsAction = numpy.dot(self.valsMatrix, self.fieldIn)[:]
    return
    

  def _initData(self):
    self.data.addScalar(vtype="int", name="_numVertices",
                        value=self.numVertices,
                        format="%d")
    self.data.addScalar(vtype="int", name="_spaceDim", value=self.spaceDim,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numCells", value=self.numCells,
                        format="%d")
    self.data.addScalar(vtype="int", name="_cellDim", value=self.cellDim,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numCorners", value=
                        self.numCorners,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numQuadPts",
                        value=self.numQuadPts,
                        format="%d")
    self.data.addScalar(vtype="int", name="_fiberDim",
                        value=self.fiberDim,
                        format="%d")
    
    self.data.addArray(vtype="double", name="_vertices", values=self.vertices,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="int", name="_cells", values=self.cells,
                       format="%8d", ncols=self.numVertices)
    
    self.data.addArray(vtype="double", name="_quadPts",
                       values=self.quadPts,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_quadWts", values=self.quadWts,
                       format="%16.8e", ncols=self.numQuadPts)
    
    self.data.addArray(vtype="double", name="_basis", values=self.basis,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_basisDeriv",
                       values=self.basisDeriv,
                       format="%16.8e", ncols=self.cellDim)
    
    self.data.addArray(vtype="double", name="_fieldIn",
                       values=self.fieldIn,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="double", name="_valsAction",
                       values=self.valsAction,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="double", name="_valsMatrix",
                       values=self.valsMatrix,
                       format="%16.8e", ncols=self.spaceDim)
      
    return

  
# End of file 
