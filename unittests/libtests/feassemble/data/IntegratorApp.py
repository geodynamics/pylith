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
    ## @li \b data Manager for output data.
    ## @li \b mesh Mesh information.
    ## @li \b quadrature Quadrature information.
    ## @li \b material Material information.
    ## @li \b solution Solution information.

    import pyre.inventory

    from pylith.utils.CppData import CppData
    data = pyre.inventory.facility("data", factory=CppData)
    data.meta['tip'] = "Data output manager."

    mesh = pyre.inventory.facility("mesh", family="mesh")
    mesh.meta['tip'] = "Mesh information."

    quadrature = pyre.inventory.facility("quadrature", family="quadrature")
    quadrature.meta['tip'] = "Quadrature information."

    material = pyre.inventory.facility("material", family="material")
    material.meta['tip'] = "Material information."

    solution = pyre.inventory.facility("solution", family="solution")
    solution.meta['tip'] = "Solution information."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="integratorapp"):
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
    self.verticesRef = None

    # This quadrature information is set by quadrature.calculateBasis()
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
    self._collectData()
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
    self.mesh = self.inventory.mesh
    self.quadrature = self.inventory.quadrature
    self.material = self.inventory.material
    self.solution = self.inventory.solution
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
    self.quadPts = self.quadrature.quadPtsRef
    self.quadWts = self.quadrature.quadWts
    (self.basis, self.basisDeriv) = self.quadrature.calculateBasis()

    # Material information
    self.matType = self.material.type
    self.matDBFilename = self.material.dbFilename
    self.matId = self.material.id
    self.matLabel = self.material.label
    self.density = self.material.density
    self.lameMu = self.material.lameMu
    self.lameLambda = self.material.lameLambda

    # Solution information
    self.dt = self.solution.dt
    self.fieldTpdt = self.solution.fieldTpdt
    self.fieldT = self.solution.fieldT
    self.fieldTmdt = self.solution.fieldTmdt
    return
  

  def _calculateResidual(self):
    """
    Calculate contribution to residual of operator for integrator.
    """
    raise NotImplementedError("Not implemented in abstract base class.")
    return


  def _calculateJacobian(self):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.
    """
    raise NotImplementedError("Not implemented in abstract base class.")
    return
    

  def _initData(self):
    # Mesh information
    self.data.addScalar(vtype="int", name="_spaceDim", value=self.spaceDim,
                        format="%d")
    self.data.addScalar(vtype="int", name="_cellDim", value=self.cellDim,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numVertices",
                        value=self.numVertices,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numCells", value=self.numCells,
                        format="%d")
    self.data.addArray(vtype="double", name="_vertices", values=self.vertices,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="int", name="_cells", values=self.cells,
                       format="%d", ncols=self.numBasis)    
    self.data.addArray(vtype="double", name="_verticesRef", values=self.verticesRef,
                       format="%16.8e", ncols=self.cellDim)

    # Quadrature information
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
                        format='"%s"')
    self.data.addScalar(vtype="char*", name="_matDBFilename",
                        value=self.matDBFilename, format='"%s"')
    self.data.addScalar(vtype="int", name="_matId",
                        value=self.matId, format="%d")
    self.data.addScalar(vtype="char*", name="_matLabel",
                        value=self.matLabel, format='"%s"')

    # Input files
    self.data.addScalar(vtype="double", name="_dt", value=self.dt,
                        format="%16.8e")
    self.data.addArray(vtype="double", name="_fieldTpdt",
                       values=self.fieldTpdt,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="double", name="_fieldT",
                       values=self.fieldT,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="double", name="_fieldTmdt",
                       values=self.fieldTmdt,
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
