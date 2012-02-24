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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pyre/meshio/DataWriterHDF5ExtSubSubMesh.py
##
## @brief Python object for writing finite-element data to HDF5Ext file.

from DataWriterHDF5Ext import DataWriterHDF5Ext
from meshio import SubSubMeshDataWriterHDF5Ext as ModuleDataWriterHDF5Ext

# DataWriterHDF5ExtSubSubMesh class
class DataWriterHDF5ExtSubSubMesh(DataWriterHDF5Ext, ModuleDataWriterHDF5Ext):
  """
  Python object for writing finite-element data to HDF5Ext file.

  Inventory

  Factory: output_data_writer
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="datawriterhdf5submesh"):
    """
    Constructor.
    """
    DataWriterHDF5Ext.__init__(self, name)
    ModuleDataWriterHDF5Ext.__init__(self)
    return


  def initialize(self, normalizer):
    """
    Initialize writer.
    """
    DataWriterHDF5Ext.initialize(self, normalizer)

    timeScale = normalizer.timeScale()

    ModuleDataWriterHDF5Ext.filename(self, self.filename)
    ModuleDataWriterHDF5Ext.timeScale(self, timeScale.value)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def data_writer():
  """
  Factory associated with DataWriter.
  """
  return DataWriterHDF5ExtSubSubMesh()


# End of file 
