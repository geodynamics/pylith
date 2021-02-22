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

## @file postproc/princaxes

## @brief Python application to compute principal axes for a tensor respresented
## as a vector. Information is read from a VTK file and a VTK file with the
## same dimensions is output.

import math
import numpy
import os
import re
import glob
from pythia.pyre.units.time import s

from pythia.pyre.applications.Script import Script as Application

class PrincAxes(Application):
  """
  Python application to compute principal axes for a tensor respresented
  as a vector. Information is read from a VTK file and a VTK file with the
  same dimensions is output.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing PrincAxes facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing PrincAxes facilities and properties.
    ##
    ## \b Properties
    ## @li \b vtk_input_file  Name of VTK input file.
    ## @li \b vtk_output_file Name of VTK output file.
    ## @li \b vtk_tensor_index Index of desired VTK field array.
    ## @li \b vtk_tensor_components_order Indices of xx,yy,zz,xy,yz,xz.
    ## @li \b add_regional_field Add regional stress/strain field?
    ## @li \b regional_sigma1 Value and direction cosines for sigma_1.
    ## @li \b regional_sigma2 Value and direction cosines for sigma_2.
    ## @li \b regional_sigma3 Value and direction cosines for sigma_3.

    import pythia.pyre.inventory

    vtkInputFile = pythia.pyre.inventory.str("vtk_input_file",
                                          default="stress_t0001.vtk")
    vtkInputFile.meta['tip'] = "Name of VTK input file."

    vtkOutputFile = pythia.pyre.inventory.str("vtk_output_file", default="output.vtk")
    vtkOutputFile.meta['tip'] = "Name of VTK output file."

    vtkTensorIndex = pythia.pyre.inventory.int("vtk_tensor_index", default=1)
    vtkTensorIndex.meta['tip'] = "Index of desired VTK field array."

    vtkTensorComponentsOrder = pythia.pyre.inventory.list(
      "vtk_tensor_components_order", default=[0, 1, 2, 3, 4, 5])
    vtkTensorComponentsOrder.meta['tip'] = "Indices of xx, yy, zz, xy, yz, xz."

    addRegionalField = pythia.pyre.inventory.bool("add_regional_field", default=False)
    addRegionalField.meta['tip'] = "Add regional field?"

    regionalSigma1 = pythia.pyre.inventory.list("regional_sigma1",
                                         default=[-1.5e7, 1.0, 0.0, 0.0])
    regionalSigma1.meta['tip'] = "Value and direction cosines of sigma1."

    regionalSigma2 = pythia.pyre.inventory.list("regional_sigma2",
                                         default=[-1.0e7, 0.0, 0.0, 1.0])
    regionalSigma2.meta['tip'] = "Value and direction cosines of sigma2."

    regionalSigma3 = pythia.pyre.inventory.list("regional_sigma3",
                                         default=[-5.0e6, 0.0, 1.0, 0.0])
    regionalSigma3.meta['tip'] = "Value and direction cosines of sigma3."
    
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="princaxes"):
    Application.__init__(self, name)

    self.cells = None
    self.vertArray = None
    self.spaceDim = 0
    self.cellType = ""

    self.numTensorPoints = 0
    self.tensorSorted = None

    self.minPrincAxis = None
    self.intPrincAxis = None
    self.maxPrincAxis = None
    self.minEigenvalue = None
    self.intEigenvalue = None
    self.maxEigenvalue = None

    self.regionalField = numpy.zeros((3, 3), dtype=numpy.float64)
    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readVtkFile()
    if (self.addRegionalField):
      self._getRegionalField()
    self._getPrincAxes()
    self._writeVtkFile()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    # import pdb
    # pdb.set_trace()

    # File info.
    self.vtkInputFile = self.inventory.vtkInputFile
    self.vtkOutputFile = self.inventory.vtkOutputFile

    # Index information
    self.vtkTensorIndex = self.inventory.vtkTensorIndex
    self.vtkTensorComponentsOrder = self.inventory.vtkTensorComponentsOrder

    # Regional field
    s1 = float(self.inventory.regionalSigma1[0])
    s2 = float(self.inventory.regionalSigma1[0])
    s3 = float(self.inventory.regionalSigma1[0])
    sVec = numpy.array([s1, s2, s3], dtype=numpy.float64)
    self.regionalSigma = numpy.diag(sVec)
    self.regionalAxes = numpy.zeros((3, 3), dtype=numpy.float64)
    for i in range(3):
      self.regionalAxes[0, i] = float(self.inventory.regionalSigma1[i+1])
      self.regionalAxes[1, i] = float(self.inventory.regionalSigma2[i+1])
      self.regionalAxes[2, i] = float(self.inventory.regionalSigma3[i+1])

    return


  def _getRegionalField(self):
    """
    Function to transform regional field from principal axes to mesh
    coordinates.
    """
    t1 = numpy.dot(self.regionalAxes, self.regionalSigma)
    self.regionalField = numpy.dot(t1, numpy.transpose(self.regionalAxes))
    return
      

  def _readVtkFile(self):
    """
    Function to read tensor from a file and store the info in a numpy array.
    """
    from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
    from enthought.tvtk.api import tvtk

    reader = VTKFileReader()
    reader.initialize(self.vtkInputFile)
    data = reader.outputs[0]

    # Get vertex and cell info.
    cellVtk = data.get_cells()
    numCells = cellVtk.number_of_cells
    cellArray = cellVtk.to_array()
    self.cells = tvtk.CellArray()
    self.cells.set_cells(numCells, cellArray)
    self.vertArray = data._get_points().to_array()
    self.cellType = data.get_cell_type(0)
    (numVerts, self.spaceDim) = self.vertArray.shape

    # Get cell fields and extract tensor.
    cellData = data._get_cell_data()
    numCellDataArrays = cellData._get_number_of_arrays()
    tensor = cellData.get_array(self.vtkTensorIndex).to_array()
    (self.numTensorPoints, numCols) = tensor.shape
    
    sxx = tensor[:, self.vtkTensorComponentsOrder[0]]
    syy = tensor[:, self.vtkTensorComponentsOrder[1]]
    szz = tensor[:, self.vtkTensorComponentsOrder[2]]
    sxy = tensor[:, self.vtkTensorComponentsOrder[3]]
    syz = tensor[:, self.vtkTensorComponentsOrder[4]]
    sxz = tensor[:, self.vtkTensorComponentsOrder[5]]
    self.tensorSorted = numpy.column_stack((sxx, syy, szz, sxy, syz, sxz))

    return


  def _getPrincAxes(self):
    """
    Function to loop over integration points and compute principal axes for
    each point.
    """
    # Create empty arrays for each principal axis and eigenvalue.
    self.minPrincAxis = numpy.empty((self.numTensorPoints, self.spaceDim),
                                    dtype=numpy.float64)
    self.intPrincAxis = numpy.empty((self.numTensorPoints, self.spaceDim),
                                    dtype=numpy.float64)
    self.maxPrincAxis = numpy.empty((self.numTensorPoints, self.spaceDim),
                                    dtype=numpy.float64)
    self.minEigenValue = numpy.empty(self.numTensorPoints, dtype=numpy.float64)
    self.intEigenValue = numpy.empty(self.numTensorPoints, dtype=numpy.float64)
    self.maxEigenValue = numpy.empty(self.numTensorPoints, dtype=numpy.float64)
    # Loop over integration points.
    for point in range(self.numTensorPoints):
      tensor = self.tensorSorted[point,:]
      tensorOrdered, eigenValuesOrdered = self._compPrincAxes(tensor)
      self.minPrincAxis[point,:] = tensorOrdered[:, 0]
      self.intPrincAxis[point,:] = tensorOrdered[:, 1]
      self.maxPrincAxis[point,:] = tensorOrdered[:, 2]
      self.minEigenValue[point] = eigenValuesOrdered[0]
      self.intEigenValue[point] = eigenValuesOrdered[1]
      self.maxEigenValue[point] = eigenValuesOrdered[2]

    return
  

  def _compPrincAxes(self, tensor):
    """
    Function to compute 3D principal axes and sort them.
    """
    tensorMat = numpy.array([(tensor[0], tensor[3], tensor[5]),
                             (tensor[3], tensor[1], tensor[4]),
                             (tensor[5], tensor[4], tensor[2])],
                            dtype=numpy.float64)
    tensorMat += self.regionalField
    (eigenValue, princAxes) = numpy.linalg.eigh(tensorMat)
    idx = eigenValue.argsort()
    eigenValuesOrdered = eigenValue[idx]
    princAxesOrdered = princAxes[:, idx]
    # tensorOrdered = numpy.empty_like(princAxesOrdered)
    # tensorOrdered[0,:] = eigenValuesOrdered[0] * princAxesOrdered[0,:]
    # tensorOrdered[1,:] = eigenValuesOrdered[1] * princAxesOrdered[1,:]
    # tensorOrdered[2,:] = eigenValuesOrdered[2] * princAxesOrdered[2,:]
    return princAxesOrdered, eigenValuesOrdered
  

  def _writeVtkFile(self):
    """
    Function to write out vertex and cell info along with principal axes
    computed as vectors.
    """
    from enthought.tvtk.api import tvtk

    # Set up mesh info for VTK file.
    mesh = tvtk.UnstructuredGrid(points=self.vertArray)
    mesh.set_cells(self.cellType, self.cells)

    # Add scalar fields.
    minEigenName = "min_eigenvalue"
    intEigenName = "int_eigenvalue"
    maxEigenName = "max_eigenvalue"
    mesh.cell_data.scalars = self.minEigenValue
    mesh.cell_data.scalars.name = minEigenName
    s2 = mesh.cell_data.add_array(self.intEigenValue)
    mesh.cell_data.get_array(s2).name = intEigenName
    s3 = mesh.cell_data.add_array(self.maxEigenValue)
    mesh.cell_data.get_array(s3).name = maxEigenName
    mesh.update()

    # Add vector fields and write VTK file
    minAxisName = "min_principal_axis"
    intAxisName = "int_principal_axis"
    maxAxisName = "max_principal_axis"
    mesh.cell_data.vectors = self.minPrincAxis
    mesh.cell_data.vectors.name = minAxisName
    v2 = mesh.cell_data.add_array(self.intPrincAxis)
    mesh.cell_data.get_array(v2).name = intAxisName
    v3 = mesh.cell_data.add_array(self.maxPrincAxis)
    mesh.cell_data.get_array(v3).name = maxAxisName
    mesh.update()
    w = tvtk.UnstructuredGridWriter(file_name=self.vtkOutputFile,
		    input=mesh)
    w.write()

    return
  

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = PrincAxes()
  app.run()

# End of file
