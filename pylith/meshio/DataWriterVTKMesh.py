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
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pyre/meshio/DataWriterVTKMesh.py
##
## @brief Python object for writing finite-element data to VTK file.

from DataWriterVTK import DataWriterVTK
from meshio import MeshDataWriterVTK as ModuleDataWriterVTK

# DataWriterVTKMesh class
class DataWriterVTKMesh(DataWriterVTK, ModuleDataWriterVTK):
  """
  Python object for writing finite-element data to VTK file.

  Inventory

  Factory: output_data_writer
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="datawritervtkmesh"):
    """
    Constructor.
    """
    DataWriterVTK.__init__(self, name)
    ModuleDataWriterVTK.__init__(self)
    return


  def initialize(self, normalizer):
    """
    Initialize writer.
    """
    DataWriterVTK.initialize(self, normalizer)
    
    timeScale = normalizer.timeScale()
    timeConstantN = normalizer.nondimensionalize(self.timeConstant,
                                                 timeScale)

    ModuleDataWriterVTK.filename(self, self.filename)
    ModuleDataWriterVTK.timeScale(self, timeScale.value)
    ModuleDataWriterVTK.timeFormat(self, self.timeFormat)
    ModuleDataWriterVTK.timeConstant(self, timeConstantN)
    ModuleDataWriterVTK.precision(self, self.precision)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def data_writer():
  """
  Factory associated with DataWriter.
  """
  return DataWriterVTKMesh()


# End of file 
