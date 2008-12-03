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

## @file pyre/meshio/DataWriterVTK.py
##
## @brief Python object for writing finite-element data to VTK file.
##
## Factory: output_data_writer

from DataWriter import DataWriter

# DataWriterVTK class
class DataWriterVTK(DataWriter):
  """
  Python object for writing finite-element data to VTK file.

  Factory: output_data_writer
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(DataWriter.Inventory):
    """
    Python object for managing DataWriterVTK facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing DataWriterVTK facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Name of VTK file.
    ## @li \b time_format C style format string for time stamp in filename.
    ## @li \b time_constant Value used to normalize time stamp in filename.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="output.vtk")
    filename.meta['tip'] = "Name of VTK file."

    timeFormat = pyre.inventory.str("time_format", default="%f")
    timeFormat.meta['tip'] = "C style format string for time stamp in filename."

    from pyre.units.time import second
    timeConstant = pyre.inventory.dimensional("time_constant",
                                              default=1.0*second,
                              validator=pyre.inventory.greater(0.0*second))
    timeConstant.meta['tip'] = "Values used to normalize time stamp in filename."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solutioniovtk"):
    """
    Constructor.
    """
    DataWriter.__init__(self, name)
    return


  def initialize(self, normalizer):
    """
    Initialize writer.
    """
    DataWriter.initialize(self, normalizer)

    # Nondimensionalize
    timeScale = normalizer.timeScale()
    self.timeConstant = normalizer.nondimensionalize(self.timeConstant,
                                                     timeScale)
    
    self.cppHandle.filename = self.filename
    self.cppHandle.timeFormat = self.timeFormat
    self.cppHandle.timeConstant = self.timeConstant
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    DataWriter._configure(self)
    self.filename = self.inventory.filename
    self.timeFormat = self.inventory.timeFormat
    self.timeConstant = self.inventory.timeConstant
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.meshio.meshio as bindings
      self.cppHandle = bindings.DataWriterVTK()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def output_data_writer():
  """
  Factory associated with DataWriterVTK.
  """
  return DataWriterVTK()


# End of file 
