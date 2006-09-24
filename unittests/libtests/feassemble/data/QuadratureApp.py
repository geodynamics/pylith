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

## @file unittests/libtests/feassemble/data/QuadratureApp.py

## @brief Python application for generating C++ data files for testing
## C++ quadrature objects.

from pyre.applications.Script import Script

# QuadratureApp class
class QuadratureApp(Script):
  """
  Python application for generating C++ data files for testing C++
  quadrature objects.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Script.Inventory):
    """Python object for managing PyLithApp facilities and properties."""

    ## @class Inventory
    ## Python object for managing PyLithApp facilities and properties.
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

  def __init__(self, name="quadratureapp"):
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

    self.quadPtsRef = None
    self.quadWts = None
    self.vertices = None
    self.cells = None

    self.basis = None
    self.basisDeriv = None
    self.quadPts = None
    self.jacobian = None
    self.jacobianDet = None
    self.jacobianInv = None
    return


  def main(self):
    """
    Run the application.
    """
    self._calculate()
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
    
    self.data.addArray(vtype="double", name="_vertices", values=self.vertices,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="int", name="_cells", values=self.cells,
                       format="%8d", ncols=self.numVertices)
    
    self.data.addArray(vtype="double", name="_quadPtsRef",
                       values=self.quadPtsRef,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_quadWts", values=self.quadWts,
                       format="%16.8e", ncols=self.numQuadPts)
    
    self.data.addArray(vtype="double", name="_basis", values=self.basis,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_basisDeriv",
                       values=self.basisDeriv,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_quadPts",
                       values=self.quadPts,
                       format="%16.8e", ncols=self.spaceDim)
    
    self.data.addArray(vtype="double", name="_jacobian",
                       values=self.jacobian,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="double", name="_jacobianDet",
                       values=self.jacobianDet,
                       format="%16.8e", ncols=self.numQuadPts)
    self.data.addArray(vtype="double", name="_jacobianInv",
                       values=self.jacobianInv,
                       format="%16.8e", ncols=self.cellDim)
      
    return

  
  def _calculate(self):
    """
    Calculate basis functions, derivatives, and Jacobian information
    at quadrature points.
    """
    raise NotImplementedError


# End of file 
