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

from DataWriter import DataWriter

# DataWriterVTK class
class DataWriterVTK(DataWriter):
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

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="datawritervtk"):
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
    return


# End of file 
