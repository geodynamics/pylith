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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pyre/meshio/DataWriterHDF5Mesh.py
##
## @brief Python object for writing finite-element data to HDF5 file.

from DataWriterHDF5 import DataWriterHDF5
from meshio import MeshDataWriterHDF5 as ModuleDataWriterHDF5

# DataWriterHDF5Mesh class
class DataWriterHDF5Mesh(DataWriterHDF5, ModuleDataWriterHDF5):
  """
  Python object for writing finite-element data to HDF5 file.

  Inventory

  Factory: output_data_writer
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="datawriterhdf5mesh"):
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
    
    ModuleDataWriterHDF5.filename(self, self.filename)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def data_writer():
  """
  Factory associated with DataWriter.
  """
  return DataWriterHDF5Mesh()


# End of file 
