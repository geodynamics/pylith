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

## @file pyre/meshio/DataWriterHDF5SubMesh.py
##
## @brief Python object for writing finite-element data to HDF5 file.

from DataWriterHDF5 import DataWriterHDF5
from meshio import SubMeshDataWriterHDF5 as ModuleDataWriterHDF5

# DataWriterHDF5SubMesh class
class DataWriterHDF5SubMesh(DataWriterHDF5, ModuleDataWriterHDF5):
  """
  Python object for writing finite-element data to HDF5 file.

  Inventory

  Factory: output_data_writer
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="datawriterhdf5submesh"):
    """
    Constructor.
    """
    DataWriterHDF5.__init__(self, name)
    ModuleDataWriterHDF5.__init__(self)
    return


  def initialize(self, normalizer):
    """
    Initialize writer.
    """
    DataWriterHDF5.initialize(self, normalizer)

    timeScale = normalizer.timeScale()
    
    ModuleDataWriterHDF5.filename(self, self.filename)
    ModuleDataWriterHDF5.timeScale(self, timeScale.value)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def data_writer():
  """
  Factory associated with DataWriter.
  """
  return DataWriterHDF5SubMesh()

 
# End of file 
