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

import numpy

# QuadratureApp class
class QuadratureApp(Script):
  """
  Python application for generating C++ data files for testing C++
  quadrature objects.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Script.Inventory):
    """Python object for managing QuadratureApp facilities and properties."""

    ## @class Inventory
    ## Python object for managing QuadratureApp facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b data Manager for output data.
    ## @li \b mesh Mesh information.
    ## @li \b quadrature Quadrature information.

    import pyre.inventory

    from pylith.utils.CppData import CppData
    data = pyre.inventory.facility("data", factory=CppData)
    data.meta['tip'] = "Manager for output data."

    mesh = pyre.inventory.facility("mesh", family="mesh")
    mesh.meta['tip'] = "Mesh information."

    quadrature = pyre.inventory.facility("quadrature", family="quadrature")
    quadrature.meta['tip'] = "Quadrature information."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadratureapp"):
    """
    Constructor.
    """
    Script.__init__(self, name)

    # Mesh information
    self.spaceDim = None
    self.cellDim = None
    self.numVertices = None
    self.numCells = None
    self.vertices = None
    self.cells = None

    # Reference cell information
    self.numBasis = None
    self.numQuadPts = None
    self.verticesRef = None
    self.quadPtsRef = None
    self.quadWts = None
    self.basis = None
    self.basisDeriv = None

    # Computed quadrature information
    self.quadPts = None
    self.jacobian = None
    self.jacobianDet = None
    self.jacobianInv = None
    return


  def main(self):
    """
    Run the application.
    """
    self._collectData()

    (self.basis, self.basisDeriv) = self.quadrature.calculateBasis()
    self.quadPts = numpy.dot(self.basis, self.vertices)

    import feutils
    (self.jacobian, self.jacobianInv, self.jacobianDet) = \
                    feutils.calculateJacobian(self.quadrature, self.vertices)

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
    self.mesh = self.inventory.mesh
    self.quadrature = self.inventory.quadrature
    return


  def _collectData(self):
    """
    Collect data we need from data objects.
    """
    # Mesh information
    self.spaceDim = self.mesh.spaceDim
    self.cellDim = self.mesh.cellDim
    self.numVertices = self.mesh.numVertices
    self.numCells = self.mesh.numCells
    self.vertices = self.mesh.vertices
    self.cells = self.mesh.cells
    self.verticesRef = self.mesh.verticesRef

    # Quadrature information
    self.numBasis = self.quadrature.numBasis
    self.numQuadPts = self.quadrature.numQuadPts
    self.quadPtsRef = self.quadrature.quadPtsRef
    self.quadWts = self.quadrature.quadWts
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
    self.data.addScalar(vtype="int", name="_numBasis", value=
                        self.numBasis,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numQuadPts",
                        value=self.numQuadPts,
                        format="%d")
    
    self.data.addArray(vtype="double", name="_vertices", values=self.vertices,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="int", name="_cells", values=self.cells,
                       format="%8d", ncols=self.numVertices)
    
    self.data.addArray(vtype="double", name="_verticesRef",
                       values=self.verticesRef,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_quadPtsRef",
                       values=self.quadPtsRef,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_quadWts", values=self.quadWts,
                       format="%16.8e", ncols=self.numQuadPts)
    self.data.addArray(vtype="double", name="_quadPts",
                       values=self.quadPts,
                       format="%16.8e", ncols=self.spaceDim)
        
    self.data.addArray(vtype="double", name="_basis",
                       values=self.basis,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_basisDeriv",
                       values=self.basisDeriv,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_jacobian",
                       values=self.jacobian,
                       format="%16.8e", ncols=self.cellDim)
    self.data.addArray(vtype="double", name="_jacobianDet",
                       values=self.jacobianDet,
                       format="%16.8e", ncols=self.numQuadPts)
    self.data.addArray(vtype="double", name="_jacobianInv",
                       values=self.jacobianInv,
                       format="%16.8e", ncols=self.spaceDim)
      
    return

  
# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = QuadratureApp()
  app.run()


# End of file 
