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

## @file greensfns/gfgen

## @brief Python application to set up impulses for generating Green's functions
## using PyLith. This application creates the necessary spatialdb files, a
## metadata file describing each impulse, and a PyLith .cfg file to set up the
## impulses.

import math
import numpy
import os
import re
import glob
from pyre.units.time import s

from pyre.applications.Script import Script as Application

class GfGen(Application):
  """
  Python application to set up impulses for generating Green's functions
  using PyLith. This application creates the necessary spatialdb files, a
  metadata file describing each impulse, and a PyLith .cfg file to set up the
  impulses.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing GfGen facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing GfGen facilities and properties.
    ##
    ## \b Properties
    ## @li \b fault_info_file  Name of VTK file containing fault info.
    ## @li \b spatialdb_output_root Root name for output spatialdb files.
    ## @li \b metadata_output_file Name of output file containing metadata.
    ## @li \b response_output_root Root name for files containing responses.
    ## @li \b impulse_type Type of impulse to be applied.
    ## @li \b impulse_value Amount of impulse to apply.
    ## @li \b timestamp_width Width of timestamp field.
    ##
    ## \b Facilities
    ## @li \b geometry  Geometry for output database.
    ## @li \b iohandler  Object for writing database.

    import pyre.inventory

    faultInfoFile = pyre.inventory.str("fault_info_file",
                                       default="fault_info.vtk")
    faultInfoFile.meta['tip'] = "Name of VTK file containing fault info."

    spatialdbOutputRoot = pyre.inventory.str("spatialdb_output_root",
                                             default="impulse.spatialdb")
    spatialdbOutputRoot.meta['tip'] = "Root name for output spatialdb files."

    metadataOutputFile = pyre.inventory.str("metadata_output_file",
                                       default="impulse_description.txt")
    metadataOutputFile.meta['tip'] = "Name of output file containing metadata."

    responseOutputRoot = pyre.inventory.str("response_output_root",
                                            default="response.vtk")
    responseOutputRoot.meta['tip'] = "Root name for files containing responses."

    impulseType = pyre.inventory.str("impulse_type",
                                     default="left-lateral-slip",
                                     validator=pyre.inventory.choice([
      "left-lateral-slip","reverse-slip","fault-opening"]))
    impulseType.meta['tip'] = "Type of impulse to apply."

    impulseValue = pyre.inventory.dimensional("impulse_value", default=1.0*m)
    impulseValue.meta['tip'] = "Impulse value."

    timestampWidth = pyre.inventory.int("timestamp_width", default=4)
    timestampWidth.meta['tip'] = "Width of timestamp field in spatialdb files."

    from spatialdata.spatialdb.generator.Geometry import Geometry
    geometry = pyre.inventory.facility("geometry", family="geometry",
                                       factory=Geometry)
    geometry.meta['tip'] = "Geometry for output database."

    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    iohandler = pyre.inventory.facility("iohandler", family="simpledb_io",
                                        factory=SimplIOAscii)
    iohandler.meta['tip'] = "Object for writing database."
    
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="gfgen"):
    Application.__init__(self, name)

    self.numFaultVertices = 0
    self.numFaultCells = 0
    self.spaceDim = 0
    self.cellType = ""

    self.faultCoords = None
    self.faultNormals = None
    self.integratedSlip = None

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readFaultInfo()
    self._makeSpatialdb()
    self._makeConfig()
    self._makeMetadata()
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
    self.faultInfoFile = self.inventory.faultInfoFile
    self.spatialdbOutputRoot = self.inventory.spatialdbOutputRoot
    self.metadataOutputFile = self.inventory.metadataOutputFile
    self.responseOutputRoot = self.inventory.responseOutputRoot

    # Impulse information
    self.impulseType = self.inventory.impulseType
    self.impulseValue = self.inventory.impulseValue.value
    self.timestampWidth = self.inventory.timestampWidth

    # Spatialdb output facilities
    self.geometry = self.inventory.geometry
    self.iohandler = self.inventory.iohandler

    return
      

  def _readFaultInfo(self):
    """
    Function to read fault information from VTK file.
    """
    from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
    from enthought.tvtk.api import tvtk

    reader = VTKFileReader()
    reader.initialize(self.faultInfoFile)
    data = reader.outputs[0]

    # Get cell and vertex info.
    cellVtk = data.get_cells()
    self.numFaultCells = cellVtk.number_of_cells
    self.faultCellArray = cellVtk.to_array()
    self.faultVertices = data._get_points().to_array()
    self.cellType = data.get_cell_type(0)
    (self.numFaultVertices, self.spaceDim) = self.faultVertices.shape

    # Get vertex fields and extract normal vectors.
    vertData = data._get_point_data()
    numVertDataArrays = vertData._get_number_of_arrays()
    for vertDataArray in range(numVertDataArrays):
      arrayName = vertData.get_array_name(vertDataArray)
      if (arrayName == "normal_dir"):
        self.faultNormals = vertData.get_array(vertDataArray).to_array()
        break

    return


  def _makeSpatialdb(self):
    """
    Function to generate a set of spatial databases (one for each impulse).
    """

    # Create empty arrays for impulse values.
    # Only array1 will be modified.
    array1 = numpy.zeros( (self.numFaultVertices,), dtype=numpy.float64)
    array2 = numpy.zeros( (self.numFaultVertices,), dtype=numpy.float64)

    # Set up info for arrays that won't be modified.
    if (self.impulseType == "left-lateral-slip"):
      info2 = {'name': "fault-opening",
               'units': "m",
               'data': array2.flatten()}
      if (self.spaceDim == 3):
        info3 = {'name': "reverse-slip",
                 'units': "m",
                 'data': array2.flatten()}
    elif (self.impulseType == "fault-opening"):
      info2 = {'name': "left-lateral-slip",
               'units': "m",
               'data': array2.flatten()}
      if (self.spaceDim == 3):
        info3 = {'name': "reverse-slip",
                 'units': "m",
                 'data': array2.flatten()}
    elif (self.impulseType == "reverse-slip"):
      info2 = {'name': "left-lateral-slip",
               'units': "m",
               'data': array2.flatten()}
      info3 = {'name': "reverse-slip",
               'units': "m",
               'data': array2.flatten()}

    # Create root output filename.
    suffIndex = self.spatialdbOutputRoot.rfind(".spatialdb")
    if (suffIndex == -1):
      outputRoot = self.spatialdbOutputRoot
    else:
      outputRoot = self.spatialdbOutputRoot[:suffIndex - 1]
    
    # Loop over impulses to generate and modify the appropriate entries.
    for impulse in range(self.numFaultVertices):

      # Set filename
      impulseNum = int(impulse)
      impulseString = repr(impulseNum).rjust(self.timestampWidth, '0')
      filename = outputRoot + "_i" + impulseString + ".spatialdb"
      self.iohandler.filename = filename

      # Modify database values.
      array1[impulse] = self.impulseValue
      if (impulse > 0):
        array1[impulse -1] = -self.impulseValue
      
      info1 = {'name': self.impulseType,
               'units': "m",
               'data': array1.flatten()}

      # Create data and write it to a database.
      if (self.spaceDim == 2):
        data = {'points': self.faultVertices,
                'coordsys': self.geometry.coordsys,
                'data_dim': self.geometry.dataDim,
                'values': [info1, info2]}
      else:
        data = {'points': self.faultVertices,
                'coordsys': self.geometry.coordsys,
                'data_dim': self.geometry.dataDim,
                'values': [info1, info2, info3]}

      self.iohandler.write(data)

    return

    
  def _getGfGen(self):
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
    for point in xrange(self.numTensorPoints):
      tensor = self.tensorSorted[point, :]
      tensorOrdered, eigenValuesOrdered = self._compGfGen(tensor)
      self.minPrincAxis[point,:] = tensorOrdered[0]
      self.intPrincAxis[point,:] = tensorOrdered[1]
      self.maxPrincAxis[point,:] = tensorOrdered[2]
      self.minEigenValue[point] = eigenValuesOrdered[0]
      self.intEigenValue[point] = eigenValuesOrdered[1]
      self.maxEigenValue[point] = eigenValuesOrdered[2]

    return
  

  def _compGfGen(self, tensor):
    """
    Function to compute 3D principal axes, sort them, and multiply by
    corresponding eigenvalue.
    """
    tensorMat = numpy.array([(tensor[0], tensor[3], tensor[5]),
                             (tensor[3], tensor[1], tensor[4]),
                             (tensor[5], tensor[4], tensor[2])],
                            dtype=numpy.float64)
    (eigenValue, princAxes) = numpy.linalg.eigh(tensorMat)
    idx = eigenValue.argsort()
    eigenValuesOrdered = eigenValue[idx]
    princAxesOrdered = princAxes[:,idx]
    tensorOrdered = numpy.empty_like(princAxesOrdered)
    tensorOrdered[0,:] = eigenValuesOrdered[0] * princAxesOrdered[0,:]
    tensorOrdered[1,:] = eigenValuesOrdered[1] * princAxesOrdered[1,:]
    tensorOrdered[2,:] = eigenValuesOrdered[2] * princAxesOrdered[2,:]
    return tensorOrdered, eigenValuesOrdered
  

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
  app = GfGen()
  app.run()

# End of file
