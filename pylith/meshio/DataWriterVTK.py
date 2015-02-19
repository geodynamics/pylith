#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pyre/meshio/DataWriterVTK.py
##
## @brief Python object for writing finite-element data to VTK file.

from DataWriter import DataWriter
from meshio import DataWriterVTK as ModuleDataWriterVTK

# DataWriterVTK class
class DataWriterVTK(DataWriter, ModuleDataWriterVTK):
  """
  Python object for writing finite-element data to VTK file.

  Inventory

  \b Properties
  @li \b filename Name of VTK file.
  @li \b time_format C style format string for time stamp in filename.
  @li \b time_constant Value used to normalize time stamp in filename.
  
  \b Facilities
  @li None
  """

  # INVENTORY //////////////////////////////////////////////////////////

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

  precision = pyre.inventory.int("float_precision", default=6,
                                 validator=pyre.inventory.greater(0))
  precision.meta['tip'] = "Precision of floating point values in output."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="datawritervtk"):
    """
    Constructor.
    """
    DataWriter.__init__(self, name)
    ModuleDataWriterVTK.__init__(self)
    return


  def initialize(self, normalizer):
    """
    Initialize writer.
    """
    DataWriter.initialize(self, normalizer, self.filename)
    
    timeScale = normalizer.timeScale()
    timeConstantN = normalizer.nondimensionalize(self.timeConstant, timeScale)

    ModuleDataWriterVTK.filename(self, self.filename)
    ModuleDataWriterVTK.timeScale(self, timeScale.value)
    ModuleDataWriterVTK.timeFormat(self, self.timeFormat)
    ModuleDataWriterVTK.timeConstant(self, timeConstantN)
    ModuleDataWriterVTK.precision(self, self.precision)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Configure object.
    """
    try:
      DataWriter._configure(self)
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring VTK output "
                       "(%s):\n%s" % (aliases, err.message))

    return


# FACTORIES ////////////////////////////////////////////////////////////

def data_writer():
  """
  Factory associated with DataWriter.
  """
  return DataWriterVTK()


# End of file 
