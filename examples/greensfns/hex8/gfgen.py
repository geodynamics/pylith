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
from pyre.units.length import m

from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
from spatialdata.geocoords.CSCart import CSCart
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
    ## @li \b slip_time_spatialdb Name of spatialdb file containing slip times.
    ## @li \b metadata_output_file Name of output file containing metadata.
    ## @li \b config_output_file Name of .cfg output file.
    ## @li \b impulse_type Type of impulse to be applied.
    ## @li \b impulse_value Amount of impulse to apply.
    ## @li \b impulse_number_width Width of impulse number field.
    ## @li \b start_impulse_num Number assigned to first impulse.
    ## @li \b exclude_file File containing vertex coordinates to exclude.
    ## @li \b distance_tol Distance tolerance to determine coincident vertices.
    ##
    ## \b Facilities
    ## @li \b geometry  Geometry for output database.

    import pyre.inventory

    faultInfoFile = pyre.inventory.str("fault_info_file",
                                       default="fault_info.vtk")
    faultInfoFile.meta['tip'] = "Name of VTK file containing fault info."

    spatialdbOutputRoot = pyre.inventory.str("spatialdb_output_root",
                                             default="impulse.spatialdb")
    spatialdbOutputRoot.meta['tip'] = "Root name for output spatialdb files."

    slipTimeSpatialdb = pyre.inventory.str("slip_time_spatialdb",
                                             default="sliptime.spatialdb")
    slipTimeSpatialdb.meta['tip'] = "Name of spatialdb file containing slip times."

    metadataOutputFile = pyre.inventory.str("metadata_output_file",
                                       default="impulse_description.txt")
    metadataOutputFile.meta['tip'] = "Name of output file containing metadata."

    configOutputFile = pyre.inventory.str("config_output_file",
                                          default="greenfn.cfg")
    configOutputFile.meta['tip'] = "Name of .cfg output file."

    impulseType = pyre.inventory.str("impulse_type",
                                     default="left-lateral-slip",
                                     validator=pyre.inventory.choice([
      "left-lateral-slip","reverse-slip","fault-opening"]))
    impulseType.meta['tip'] = "Type of impulse to apply."

    impulseValue = pyre.inventory.dimensional("impulse_value", default=1.0*m)
    impulseValue.meta['tip'] = "Impulse value."

    impulseNumberWidth = pyre.inventory.int("impulse_number_width", default=4)
    impulseNumberWidth.meta['tip'] = "Width of impulse number field."

    startImpulseNum = pyre.inventory.int("start_impulse_num", default=0)
    startImpulseNum.meta['tip'] = "Number assigned to first impulse."

    excludeFile = pyre.inventory.str("exclude_file", default="")
    excludeFile.meta['tip'] = "File containing vertex coordinates to exclude."

    distanceTol = pyre.inventory.dimensional("distance_tol", default=0.1*m)
    distanceTol.meta['tip'] = "Distance tolerance to determine coincident vertices."

    from spatialdata.spatialdb.generator.Geometry import Geometry
    geometry = pyre.inventory.facility("geometry", family="geometry",
                                       factory=Geometry)
    geometry.meta['tip'] = "Geometry for output database."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="gfgen"):
    Application.__init__(self, name)

    self.numFaultVertices = 0
    self.numExcludedVertices = 0
    self.numIncludedVertices = 0
    self.numFaultCells = 0
    self.spaceDim = 0
    self.cellType = ""

    self.excludeCoords = None
    self.numExcludeCoords = 0
    self.faultVertices = None
    self.normalDir = None
    self.strikeDir = None
    self.dipDir = None

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readExcludeInfo()
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
    self.slipTimeSpatialdb = self.inventory.slipTimeSpatialdb
    self.metadataOutputFile = self.inventory.metadataOutputFile
    self.configOutputFile = self.inventory.configOutputFile

    # Impulse information
    self.impulseType = self.inventory.impulseType
    self.impulseValue = self.inventory.impulseValue.value
    self.impulseNumberWidth = self.inventory.impulseNumberWidth
    self.startImpulseNum = self.inventory.startImpulseNum

    # Excluded vertex information
    self.excludeFile = self.inventory.excludeFile
    self.distanceTol = self.inventory.distanceTol.value

    # Spatialdb output facilities
    self.geometry = self.inventory.geometry

    return
      

  def _readExcludeInfo(self):
    """
    Function to read coordinates of vertices to exclude from impulse list.
    """
    fileExists = os.path.isfile(self.excludeFile)
    if fileExists:
      self.excludeCoords = numpy.loadtxt(self.excludeFile, dtype=numpy.float64)
      self.numExcludeCoords = self.excludeCoords.shape[0]
    else:
      print "No exclude file found!"

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
    faultCoords = data._get_points().to_array()
    self.cellType = data.get_cell_type(0)
    (self.numFaultVertices, self.spaceDim) = faultCoords.shape

    # Get vertex fields and extract normal vectors.
    vertData = data._get_point_data()
    numVertDataArrays = vertData._get_number_of_arrays()
    for vertDataArray in range(numVertDataArrays):
      arrayName = vertData.get_array_name(vertDataArray)
      if (arrayName == "normal_dir"):
        normalDir = vertData.get_array(vertDataArray).to_array()
      elif (arrayName == "strike_dir"):
        strikeDir = vertData.get_array(vertDataArray).to_array()
      elif (arrayName == "dip_dir"):
        dipDir = vertData.get_array(vertDataArray).to_array()

    # Determine which vertices to exclude, if any.
    includedVertices = []
    excludedVertices = []
    if self.numExcludeCoords == 0:
      self.numIncludedVertices = self.numFaultVertices
    else:
      for vertex in range(self.numFaultVertices):
        vertexCoords = faultCoords[vertex,:]
        vertexMatch = False
        for excludeVertex in range(self.numExcludeCoords):
          excludeCoords = self.excludeCoords[excludeVertex,:]
          diff = vertexCoords - excludeCoords
          dist =numpy.linalg.norm(diff)
          if (dist <= self.distanceTol):
            vertexMatch = True
            excludedVertices.append(vertex)
            break

        if not vertexMatch:
          includedVertices.append(vertex)
          
      self.numIncludedVertices = len(includedVertices)
      self.numExcludedVertices = len(excludedVertices)
      # Reorder arrays so excluded vertices are at the end.
      includedCoords = faultCoords[includedVertices,:]
      excludedCoords = faultCoords[excludedVertices,:]
      self.faultVertices = numpy.vstack((includedCoords, excludedCoords))
      includedNormal = normalDir[includedVertices,:]
      excludedNormal = normalDir[excludedVertices,:]
      self.normalDir = numpy.vstack((includedNormal, excludedNormal))
      includedStrike = strikeDir[includedVertices,:]
      excludedStrike = strikeDir[excludedVertices,:]
      self.strikeDir = numpy.vstack((includedStrike, excludedStrike))
      includedDip = dipDir[includedVertices,:]
      excludedDip = dipDir[excludedVertices,:]
      self.dipDir = numpy.vstack((includedDip, excludedDip))
    
    print "Number of excluded vertices:  %i" % self.numExcludedVertices
    print "Number of included vertices:  %i" % self.numIncludedVertices

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
      info3 = {'name': "fault-opening",
               'units': "m",
               'data': array2.flatten()}

    # Create root output filename.
    suffIndex = self.spatialdbOutputRoot.rfind(".spatialdb")
    outputRoot = self.spatialdbOutputRoot
    if (suffIndex != -1):
      outputRoot = self.spatialdbOutputRoot[:suffIndex - 1]

    # Set up coordinate system.
    # cs = self.geometry.coordsys()
    # cs.initialize()
    dataDim = self.spaceDim - 1
    
    # Loop over impulses to generate and modify the appropriate entries.
    impulse = self.startImpulseNum
    for vertex in range(self.numIncludedVertices):

      # Set filename
      impulseNum = int(impulse)
      impulseString = repr(impulseNum).rjust(self.impulseNumberWidth, '0')
      filename = outputRoot + "_i" + impulseString + ".spatialdb"
      writer = SimpleIOAscii()
      writer.inventory.filename = filename
      writer._configure()

      # Modify database values.
      array1[vertex] = self.impulseValue
      if (vertex > 0):
        array1[vertex - 1] = -self.impulseValue
      
      if (vertex > 1):
	array1[vertex - 2] = 0.0
      info1 = {'name': self.impulseType,
               'units': "m",
               'data': array1.flatten()}

      # Create data and write it to a database.
      if (self.spaceDim == 2):
        data = {'points': self.faultVertices,
                'coordsys': self.geometry.coordsys,
                'data_dim': dataDim,
                'values': [info1, info2]}
      else:
        data = {'points': self.faultVertices,
                'coordsys': self.geometry.coordsys,
                'data_dim': dataDim,
                'values': [info1, info2, info3]}

      writer.write(data)
      impulse += 1

    return

    
  def _makeConfig(self):
    """
    Function to create .cfg file for creating Green's functions.
    """
    f = open(self.configOutputFile, 'w')
    newLine = "\n"

    # Write top-level header
    topHeader = "# -*- Python -*-" + newLine + "[pylithapp]" + newLine + newLine
    f.write(topHeader)
    
    # Write time step information
    time = float(self.startImpulseNum) + float(self.numIncludedVertices) - 1.0
    totalTime = "total_time = " + repr(time) + "*year" + newLine
    dt  = "dt = 1.0*year" + newLine
    divider = "# ----------------------------------------------------------" + \
              newLine
    problem = divider + "# problem" + newLine + divider + \
              "[pylithapp.timedependent.implicit.time_step]" + newLine + \
              totalTime + dt
    f.write(problem)
    faults = divider + "# faults" + newLine + divider + \
             "[pylithapp.timedependent.interfaces.fault]" + newLine
    f.write(newLine)
    f.write(faults)
    eqSrcs = "eq_srcs = ["
    f.write(eqSrcs)
    comma = ","

    # Write eq_srcs list.
    impulse = self.startImpulseNum
    for vertex in range(self.numIncludedVertices):
      srcName = repr(impulse).rjust(self.impulseNumberWidth, '0')
      f.write(srcName)
      if (vertex != self.numIncludedVertices - 1):
        f.write(comma)
      impulse += 1

    srcEnd = "]" + newLine
    f.write(srcEnd)

    # Base strings for eq_src info.
    baseHeader = "[pylithapp.timedependent.interfaces.fault.eq_srcs."
    baseOrigin = "origin_time = "
    suffIndex = self.spatialdbOutputRoot.rfind(".spatialdb")
    slipRoot = self.spatialdbOutputRoot
    if (suffIndex != -1):
      slipRoot = self.spatialdbOutputRoot[:suffIndex - 1]
    baseSlip = "slip_function.slip.iohandler.filename = " + slipRoot + "_i"
    slipTime = "slip_function.slip_time.iohandler.filename = " + \
               self.slipTimeSpatialdb + newLine

    # Write info for each eq_src.
    impulse = self.startImpulseNum
    for vertex in range(self.numIncludedVertices):
      f.write(newLine)
      srcName = repr(impulse).rjust(self.impulseNumberWidth, '0')
      originTime = float(impulse)
      originString = str(originTime) + "*year\n"
      f.write(baseHeader + srcName + "]\n")
      f.write(baseOrigin + originString)
      f.write(baseSlip + srcName + ".spatialdb\n")
      f.write(slipTime)
      impulse += 1

    f.close()
    
    return
  

  def _makeMetadata(self):
    """
    Function to write out metadata file containing information for each
    impulse.
    """
    f = open(self.metadataOutputFile, 'w')
    tab = "\t"
    newLine = "\n"
    header = "Impulse #" + tab + "X-Coord" + tab + "Y-Coord" + tab + \
             "Z-Coord" + tab + "Normal-X" + tab + "Normal-Y" + tab + \
             "Normal-Z" + tab + "Strike-X" + tab + "Strike-Y" + tab + \
             "Strike-Z" + tab + "Dip-X" + tab + "Dip-Y" + tab + "Dip-Z" \
             + newLine
    f.write(header)
    impulse = self.startImpulseNum
    for vertex in range(self.numIncludedVertices):
      x = str(self.faultVertices[vertex, 0]) + tab
      y = str(self.faultVertices[vertex, 1]) + tab
      z = str(self.faultVertices[vertex, 2]) + tab
      xNorm = str(self.normalDir[vertex, 0]) + tab
      yNorm = str(self.normalDir[vertex, 1]) + tab
      zNorm = str(self.normalDir[vertex, 2]) + tab
      xStrike = str(self.strikeDir[vertex, 0]) + tab
      yStrike = str(self.strikeDir[vertex, 1]) + tab
      zStrike = str(self.strikeDir[vertex, 2]) + tab
      xDip = str(self.dipDir[vertex, 0]) + tab
      yDip = str(self.dipDir[vertex, 1]) + tab
      zDip = str(self.dipDir[vertex, 2]) + newLine
      outLine = str(impulse) + tab + x + y + z + xNorm + yNorm + zNorm + \
                xStrike + yStrike + zStrike + xDip + yDip + zDip
      f.write(outLine)
      impulse += 1

    f.close()

    return
  

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = GfGen()
  app.run()

# End of file
