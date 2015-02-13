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

## @file utils/VTKDataReader.py
##
## @brief Object for reading VTK data file.

import numpy

def has_vtk():
  if not "flag" in dir(has_vtk):
    try:
      from tvtk.api import tvtk
      has_vtk.flag = True
    except ImportError:
      print "WARNING: Cannot find Mayavi VTK interface to check output."
      print "         Tests limited to running PyLith without errors."
      print "         Install MayaVi2 " \
          "(http://code.enthought.com/projects/mayavi/)"
      print "         in order to enable verification of output."
      has_vtk.flag = False
  return has_vtk.flag
    

class VTKDataReader(object):
  """
  Object for reading VTK data file.
  """

  def __init__(self):
    from tvtk.api import tvtk
    reader = tvtk.UnstructuredGridReader()
    reader.read_all_scalars = True
    reader.read_all_vectors = True
    reader.read_all_tensors = True
    reader.read_all_fields = True
    self.reader = reader
    return


  def read(self, filename):
    """
    Read VTK file.
    """
    self.reader.file_name = filename
    self.reader.update()
    data = self.reader.get_output()
    
    data = {'vertices': self._getVertices(data),
            'cells': self._getCells(data),
            'vertex_fields': self._getVertexFields(data),
            'cell_fields': self._getCellFields(data),
            }
    return data


  def _getVertices(self, data):
    vertices = data.points.to_array()
    return vertices


  def _getCells(self, data):
    cells = data.get_cells().to_array()
    cellTypes = data.cell_types_array.to_array()
    ncells = data.number_of_cells
    cellId = cellTypes[0]
    if numpy.sum(cellTypes-cellId) != 0:
      raise ValueError("Expecting cells to all be the same type.")
    if cellId == 1: # vertex
      ncorners = 1
    elif cellId == 3: # line2
      ncorners = 2
    elif cellId == 5: # tri3
      ncorners = 3
    elif cellId == 9: # quad4
      ncorners = 4
    elif cellId == 10: # tet4
      ncorners = 4
    elif cellId == 12: # hex8
      ncorners = 8
    elif cellId == 21: # line3
      ncorners = 3
    elif cellId == 22: # tri6
      ncorners = 6
    elif cellId == 22: # quad9
      ncorners = 9
    elif cellId == 24: # tet10
      ncorners = 10
    elif cellId == 29: # hex27
      ncorners = 27
    elif cellId == 255: # unknown?
      ncorners = 1
    else:
      raise ValueError("Unknown VTK cell type '%d'." % cellId)
    cells = cells.reshape( (ncells, 1+ncorners) )[:,1:1+ncorners]
    return cells


  def _getVertexFields(self, data):
    fields = {}
    ptdata = data.point_data
    nfields = ptdata.number_of_arrays
    for i in xrange(nfields):
      name = ptdata.get_array_name(i)
      field = ptdata.get_array(name).to_array()
      fields[name] = field
    return fields


  def _getCellFields(self, data):
    fields = {}
    cdata = data.cell_data
    nfields = cdata.number_of_arrays
    for i in xrange(nfields):
      name = cdata.get_array_name(i)
      field = cdata.get_array(name).to_array()
      fields[name] = field
    return fields


# End of file 
