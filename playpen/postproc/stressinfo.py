#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

## @file postproc/stressinfo

## @brief Python application to compute several stress-related quantities
## from the stress tensor provided by PyLith. Information is read from a VTK
## file and a VTK file with the same dimensions is output.

import math
import numpy
import os
import re
import glob

from pythia.pyre.applications.Script import Script as Application

class StressInfo(Application):
  """Python application to compute several stress-related quantities
  from the stress tensor provided by PyLith. Information is read from a VTK
  file and a VTK file with the same dimensions is output.
  """
  
  class Inventory(Application.Inventory):
    """Python object for managing StressInfo facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing StressInfo facilities and properties.
    ##
    ## \b Properties
    ## @li \b vtk_input_file  Name of VTK input file.
    ## @li \b vtk_output_file Name of VTK output file.
    ## @li \b tensor_index Index of desired input VTK field array.
    ## @li \b tensor_components_order Indices of xx,yy,zz,xy,yz,xz.
    ## @li \b friction_angle Friction angle for plasticity calculation.
    ## @li \b cohesion Cohesion for plasticity calculation.

    import pythia.pyre.inventory
    from pythia.pyre.units.pressure import MPa
    from pythia.pyre.units.angle import degree

    vtkInputFile = pythia.pyre.inventory.str("vtk_input_file",
                                          default="stress_t0001.vtk")
    vtkInputFile.meta['tip'] = "Name of VTK input file."

    vtkOutputFile = pythia.pyre.inventory.str("vtk_output_file", default="output.vtk")
    vtkOutputFile.meta['tip'] = "Name of VTK output file."

    tensorIndex = pythia.pyre.inventory.int("tensor_index", default=1)
    tensorIndex.meta['tip'] = "Index of desired input VTK field array."

    tensorComponentsOrder = pythia.pyre.inventory.list("tensor_components_order",
                                                default=[0, 1, 2, 3, 4, 5])
    tensorComponentsOrder.meta['tip'] = "Indices of xx, yy, zz, xy, yz, xz."

    frictionAngle = pythia.pyre.inventory.dimensional("friction_angle",
                                               default=30.0*degree)
    frictionAngle.meta['tip'] = "Friction angle for plasticity calculation."

    cohesion = pythia.pyre.inventory.dimensional("cohesion",
                                          default=1.0*MPa)
    cohesion.meta['tip'] = "Cohesion for plasticity calculation."
    
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="stressinfo"):
    Application.__init__(self, name)

    self.alphaYield = 0.0
    self.beta = 0.0
    
    self.cells = None
    self.vertArray = None
    self.spaceDim = 0
    self.cellType = ""

    self.numTensorPoints = 0
    self.tensorSorted = None

    self.pressure = None
    self.devInvariant2 = None
    self.dpPlasPresTerm = None
    self.dpPlasStressTerm = None
    self.dpPlasYieldFunc = None
    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readVtkFile()
    self._getStressInfo()
    self._writeVtkFile()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Setup members using inventory.
    """
    Application._configure(self)
    # import pdb
    # pdb.set_trace()

    # File info.
    self.vtkInputFile = self.inventory.vtkInputFile
    self.vtkOutputFile = self.inventory.vtkOutputFile

    # Index information
    self.tensorIndex = self.inventory.tensorIndex
    self.tensorComponentsOrder = self.inventory.tensorComponentsOrder

    # Parameters
    self.frictionAngle = self.inventory.frictionAngle.value
    self.cohesion = self.inventory.cohesion.value
    sinFric = math.sin(self.frictionAngle)
    cosFric = math.cos(self.frictionAngle)
    denomFriction = math.sqrt(3.0) * (3.0 - sinFric)
    self.alphaYield = 2.0 * sinFric/denomFriction
    self.beta = 6.0 * self.cohesion * cosFric/denomFriction

    return
      

  def _readVtkFile(self):
    """Function to read tensor from a file and store the info in a numpy array.
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
    tensor = cellData.get_array(self.tensorIndex).to_array()
    (self.numTensorPoints, numCols) = tensor.shape
    
    sxx = tensor[:, self.tensorComponentsOrder[0]]
    syy = tensor[:, self.tensorComponentsOrder[1]]
    szz = tensor[:, self.tensorComponentsOrder[2]]
    sxy = tensor[:, self.tensorComponentsOrder[3]]
    syz = tensor[:, self.tensorComponentsOrder[4]]
    sxz = tensor[:, self.tensorComponentsOrder[5]]
    self.tensorSorted = numpy.column_stack((sxx, syy, szz, sxy, syz, sxz))

    return


  def _getStressInfo(self):
    """Function to loop over integration points and compute stress information for
    each point.
    """
    # Create empty arrays for each stress quantity
    self.pressure = numpy.empty(self.numTensorPoints, dtype=numpy.float64)
    self.devInvariant2 = numpy.empty(self.numTensorPoints, dtype=numpy.float64)
    self.dpPlasPresTerm = numpy.empty(self.numTensorPoints, dtype=numpy.float64)
    self.dpPlasStressTerm = numpy.empty(self.numTensorPoints,
                                        dtype=numpy.float64)
    self.dpPlasYieldFunc = numpy.empty(self.numTensorPoints,
                                       dtype=numpy.float64)
    # Loop over integration points.
    for point in range(self.numTensorPoints):
      tensor = self.tensorSorted[point,:]
      pressure, devInvariant2 = self._compStressInfo(tensor)
      self.pressure[point] = pressure
      self.devInvariant2[point] = devInvariant2
      dpPlasPresTerm = self.alphaYield * 3.0 * pressure
      self.dpPlasPresTerm[point] = dpPlasPresTerm
      dpPlasStressTerm = dpPlasPresTerm + devInvariant2
      self.dpPlasStressTerm[point] = dpPlasStressTerm
      self.dpPlasYieldFunc[point] = dpPlasStressTerm - self.beta

    return
  

  def _compStressInfo(self, tensor):
    """Function to compute the pressure and the second deviatoric invariant.
    """
    pressure = (tensor[0] + tensor[1] + tensor[2])/3.0
    dev = tensor
    dev[0] -= pressure
    dev[1] -= pressure
    dev[2] -= pressure
    scalarProd = numpy.dot(dev, dev) + \
                 dev[3] * dev[3] + dev[4] * dev[4] + dev[5] * dev[5]
    devInvariant2 = math.sqrt(0.5 * scalarProd)
    return pressure, devInvariant2
  

  def _writeVtkFile(self):
    """Function to write out vertex and cell info along with stress-related
    quantities.
    """
    from enthought.tvtk.api import tvtk

    # Set up mesh info for VTK file.
    mesh = tvtk.UnstructuredGrid(points=self.vertArray)
    mesh.set_cells(self.cellType, self.cells)

    # Add scalar fields.
    pressureName = "pressure"
    devInvariant2Name = "dev_invar2"
    dpPlasPresTermName = "dp_plas_pres_term"
    dpPlasStressTermName = "dp_plas_stress_term"
    dpPlasYieldFuncName = "dp_plas_yield_func"
    mesh.cell_data.scalars = self.pressure
    mesh.cell_data.scalars.name = pressureName
    s2 = mesh.cell_data.add_array(self.devInvariant2)
    mesh.cell_data.get_array(s2).name = devInvariant2Name
    s3 = mesh.cell_data.add_array(self.dpPlasPresTerm)
    mesh.cell_data.get_array(s3).name = dpPlasPresTermName
    s4 = mesh.cell_data.add_array(self.dpPlasStressTerm)
    mesh.cell_data.get_array(s4).name = dpPlasStressTermName
    s5 = mesh.cell_data.add_array(self.dpPlasYieldFunc)
    mesh.cell_data.get_array(s5).name = dpPlasYieldFuncName
    mesh.update()

    # Write VTK file
    w = tvtk.XMLDataSetWriter(file_name=self.vtkOutputFile, input=mesh)
    w.write()

    return
  

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = StressInfo()
  app.run()

# End of file
