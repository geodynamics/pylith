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

## @file postproc/faultinfo

## @brief Python application to compute a number of quantities on the fault
## mesh using both the fault info file and the computed fault output from
## PyLith. Information includes slip and stress vectors in global coordinates
## and the CFF change using a specified value for the effective coefficient of
## friction. Note that this application requires the output of fault tractions
## as well as all three orientation vectors.

import numpy

from pythia.pyre.applications.Script import Script as Application

class FaultInfo(Application):
  """Python application to compute a number of quantities on the fault mesh using
  both the fault info file and the computed fault output from PyLith.
  Information includes slip and stress vectors in global coordinates and the CFF
  change using a specified value for the effective coefficient of friction.
  Note that this application requires the output of fault tractions as well as
  all three orientation vectors.
  """
  
  class Inventory(Application.Inventory):
    """Python object for managing FaultInfo facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing FaultInfo facilities and properties.
    ##
    ## \b Properties
    ## @li \b fault_info_file VTK file containing fault information.
    ## @li \b fault_results_file VTK file containing computed fault results.
    ## @li \b fault_output_file VTK output file.
    ## @li \b friction_coeff Effective coefficient of friction.
    ## @li \b stress_scale_factor Stress scale factor.
    ## @li \b slip_scale_factor Fault slip scale factor.
    ## @li \b shear_direction General direction associated with positive shear.

    import pythia.pyre.inventory

    faultInfoFile = pythia.pyre.inventory.str("fault_info_file",
                                          default="fault_info.vtk")
    faultInfoFile.meta['tip'] = "VTK file containing fault information."

    faultResultsFile = pythia.pyre.inventory.str("fault_results_file",
                                          default="fault_t00000.vtk")
    faultResultsFile.meta['tip'] = "VTK file containing computed fault results."

    faultOutputFile = pythia.pyre.inventory.str("fault_output_file",
                                         default="fault_stress_t00000.vtk")
    faultOutputFile.meta['tip'] = "VTK output file."

    frictionCoeff = pythia.pyre.inventory.float("friction_coeff", default=0.2)
    frictionCoeff.meta['tip'] = "Effective coefficient of friction."

    stressScaleFactor = pythia.pyre.inventory.float("stress_scale_factor", default=1.0)
    stressScaleFactor.meta['tip'] = "Scale factor to apply to stresses."

    slipScaleFactor = pythia.pyre.inventory.float("slip_scale_factor", default=1.0)
    slipScaleFactor.meta['tip'] = "Scale factor to apply to slip."

    shearDirection = pythia.pyre.inventory.list("shear_direction",
                                         default=[1.0, 0.0, 1.0])
    shearDirection.meta['tip'] = "General direction associated with positive shear stress."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultinfo"):
    Application.__init__(self, name)

    self.numVertsPerCell = 0
    self.numCells = 0
    self.cellsArray = numpy.array([0])
    self.verticesArray = numpy.array([0])
    self.normalVec = numpy.array([0])
    self.strikeVec = numpy.array([0])
    self.dipVec = numpy.array([0])
    self.stressVec = numpy.array([0])
    self.slipVec = numpy.array([0])
    self.numVerts = 0
    self.spaceDim = 0
    self.cellType = ""
    self.readMesh = False

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readVectors()
    self._computeVectorInfo()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Setup members using inventory.
    """
    Application._configure(self)

    # Filenames
    self.faultInfoFile = self.inventory.faultInfoFile
    self.faultResultsFile = self.inventory.faultResultsFile
    self.faultOutputFile = self.inventory.faultOutputFile

    # Parameters
    self.frictionCoeff = self.inventory.frictionCoeff
    self.stressScaleFactor = self.inventory.stressScaleFactor
    self.slipScaleFactor = self.inventory.slipScaleFactor
    self.shearDirection = self.inventory.shearDirection

    return


  def _readVectors(self):
    """Get orientation vectors and stress vector.
    """
    self.normalVec = self._getVec(self.faultInfoFile, "normal_dir")
    self.strikeVec = self._getVec(self.faultInfoFile, "strike_dir")
    self.dipVec = self._getVec(self.faultInfoFile, "dip_dir")
    self.stressVec = self._getVec(self.faultResultsFile, "traction_change")
    self.slipVec = self._getVec(self.faultResultsFile, "slip")
    return


  def _getVec(self, vtkFile, vecName):
    """Function to read a vector from a file and store it in a numpy array.
    """
    from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
    from enthought.tvtk.api import tvtk

    reader = VTKFileReader()
    reader.initialize(vtkFile)
    data = reader.outputs[0]

    # Get vertex and cell info if it hasn't already been done
    if not self.readMesh:
      cellVtk = data.get_cells()
      self.numVertsPerCell = cellVtk._get_max_cell_size()
      self.numCells = cellVtk.number_of_cells
      cellArray = cellVtk.to_array()
      self.cells = tvtk.CellArray()
      self.cells.set_cells(self.numCells, cellArray)
      self.vertArray = data._get_points().to_array()
      self.cellType = data.get_cell_type(0)
      (self.numVerts, self.spaceDim) = self.vertArray.shape
      self.readMesh = True


    # Get vertex fields and extract the requested vector.
    vertData = data._get_point_data()
    numArrays = vertData._get_number_of_arrays()
    gotArray = False
    for vertDataArray in range(numArrays):
      arrayName = vertData.get_array_name(vertDataArray)
      if (arrayName == vecName):
        vector = vertData.get_array(vertDataArray).to_array()
        gotArray = True
      if gotArray:
        break
    else:
      raise IOError("Unable to find vector '%s'." % vecName)

    return vector
  

  def _computeVectorInfo(self):
    """Function to compute vectors in global coordinates and CFF.
    """
    from enthought.tvtk.api import tvtk

    # Compute stress vectors in global coordinates.
    strikeStress = numpy.tile(numpy.expand_dims(self.stressVec[:, 0], 1), (1, 3))
    strikeStressArr = self.stressScaleFactor * self.strikeVec * strikeStress
    strikeStressVec = tvtk.DoubleArray(name='left_lateral_shear')
    strikeStressVec.from_array(strikeStressArr)
    
    dipStress = numpy.tile(numpy.expand_dims(self.stressVec[:, 1], 1), (1, 3))
    dipStressArr = self.stressScaleFactor * self.dipVec * dipStress
    dipStressVec = tvtk.DoubleArray(name='up_dip_shear')
    dipStressVec.from_array(dipStressArr)

    normalStress = numpy.tile(numpy.expand_dims(self.stressVec[:, 2], 1), (1, 3))
    normalStressArr = self.stressScaleFactor * self.normalVec * normalStress
    normalStressVec = tvtk.DoubleArray(name='normal_stress')
    normalStressVec.from_array(normalStressArr)

    shearStressArr = strikeStressArr + dipStressArr
    shearStressVec = tvtk.DoubleArray(name='in_plane_shear')
    shearStressVec.from_array(shearStressArr)

    # Compute slip vectors in global coordinates.
    strikeSlip = numpy.tile(numpy.expand_dims(self.slipVec[:, 0], 1), (1, 3))
    strikeSlipArr = self.slipScaleFactor * self.strikeVec * strikeSlip
    strikeSlipVec = tvtk.DoubleArray(name='left_lateral_slip')
    strikeSlipVec.from_array(strikeSlipArr)

    dipSlip = numpy.tile(numpy.expand_dims(self.slipVec[:, 1], 1), (1, 3))
    dipSlipArr = self.slipScaleFactor * self.dipVec * dipSlip
    dipSlipVec = tvtk.DoubleArray(name='up_dip_slip')
    dipSlipVec.from_array(dipSlipArr)
    
    normalSlip = numpy.tile(numpy.expand_dims(self.slipVec[:, 2], 1), (1, 3))
    normalSlipArr = self.slipScaleFactor * self.normalVec * normalSlip
    normalSlipVec = tvtk.DoubleArray(name='normal_slip')
    normalSlipVec.from_array(normalSlipArr)

    slipArr = strikeSlipArr + dipSlipArr
    slipVec = tvtk.DoubleArray(name='in_plane_slip')
    slipVec.from_array(slipArr)

    # Compute CFF
    shearMag = self.stressScaleFactor * \
               numpy.sqrt(numpy.square(self.stressVec[:, 0]) + \
                          numpy.square(self.stressVec[:, 1]))
    
    shearDirVec = numpy.array(self.shearDirection, dtype=numpy.float64)
    shearDirMat = numpy.tile(shearDirVec, (self.numVerts, 1))
    shearSlipDot = shearStressVec * shearDirMat
    shearSlipSign = numpy.sign(numpy.sum(shearSlipDot, axis=1))
    CffArr = shearSlipSign * shearMag + \
             self.stressScaleFactor * self.frictionCoeff * self.stressVec[:, 2]
    CFF = tvtk.DoubleArray(name='CFF')
    CFF.from_array(CffArr)

    # Set up mesh info for output VTK file
    mesh = tvtk.UnstructuredGrid(points=self.vertArray)
    mesh.set_cells(self.cellType, self.cells)

    # Add computed values to mesh object.
    mesh.point_data.add_array(strikeStressVec)
    mesh.point_data.add_array(dipStressVec)
    mesh.point_data.add_array(normalStressVec)
    mesh.point_data.add_array(shearStressVec)
    mesh.point_data.add_array(strikeSlipVec)
    mesh.point_data.add_array(dipSlipVec)
    mesh.point_data.add_array(normalSlipVec)
    mesh.point_data.add_array(slipVec)
    mesh.point_data.scalars = CFF

    # Write results to VTK file
    w = tvtk.XMLDataSetWriter(file_name=self.faultOutputFile, input=mesh)
    w.write()
    
    return
    
# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = FaultInfo()
  app.run()

# End of file
