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

## @file pyre/meshio/DataWriterHDF5Ext.py
##
## @brief Python object for writing finite-element data to HDF5 file
## with datasets stored in external binary files.

from DataWriter import DataWriter
from meshio import DataWriterHDF5Ext as ModuleDataWriterHDF5Ext

# DataWriterHDF5Ext class
class DataWriterHDF5Ext(DataWriter, ModuleDataWriterHDF5Ext):
  """
  @brief Python object for writing finite-element data to HDF5 file
  with datasets stored in external binary files.

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
    ModuleDataWriterHDF5Ext.__init__(self)
    return


  def initialize(self, normalizer):
    """
    Initialize writer.
    """
    DataWriter.initialize(self, normalizer, self.filename)
    
    timeScale = normalizer.timeScale()

    ModuleDataWriterHDF5Ext.filename(self, self.filename)
    ModuleDataWriterHDF5Ext.timeScale(self, timeScale.value)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def data_writer():
  """
  Factory associated with DataWriter.
  """
  return DataWriterHDF5Ext()


# End of file 
