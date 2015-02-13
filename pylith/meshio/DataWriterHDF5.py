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

## @file pyre/meshio/DataWriterHDF5.py
##
## @brief Python object for writing finite-element data to HDF5 file.

from DataWriter import DataWriter
from meshio import DataWriterHDF5 as ModuleDataWriterHDF5

# DataWriterHDF5 class
class DataWriterHDF5(DataWriter, ModuleDataWriterHDF5):
  """
  Python object for writing finite-element data to HDF5 file.

  Inventory

  \b Properties
  @li \b filename Name of HDF5 file.
  
  \b Facilities
  @li None
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  filename = pyre.inventory.str("filename", default="output.h5")
  filename.meta['tip'] = "Name of HDF5 file."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="datawriterhdf5"):
    """
    Constructor.
    """
    DataWriter.__init__(self, name)
    ModuleDataWriterHDF5.__init__(self)
    return


  def initialize(self, normalizer):
    """
    Initialize writer.
    """
    DataWriter.initialize(self, normalizer, self.filename)

    timeScale = normalizer.timeScale()
    
    ModuleDataWriterHDF5.filename(self, self.filename)
    ModuleDataWriterHDF5.timeScale(self, timeScale.value)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def data_writer():
  """
  Factory associated with DataWriter.
  """
  return DataWriterHDF5()


# End of file 
