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

## @file utils/VTKDataReader.py
##
## @brief Object for reading VTK data file.

import numpy

def has_vtk():
  if not "flag" in dir(has_vtk):
    try:
      from enthought.tvtk.api import tvtk
      has_vtk.flag = True
    except ImportError:
      print "WARNING: Cannot find Mayavi VTK interface to check output."
      print "         Tests limited to running PyLith without errors."
      print "         Install MayaVi2 " \
          "(https://svn.enthought.com/enthought/wiki/MayaVi)"
      print "         in order to enable verification of output."
      has_vtk.flag = False
  return has_vtk.flag
    

class VTKDataReader(object):
  """
  Object for reading VTK data file.
  """

  def __init__(self):
    return


  def read(self, filename):
    from enthought.tvtk.api import tvtk
    reader = tvtk.UnstructuredGridReader()
    reader.file_name = filename
    reader.read_all_scalars = True
    reader.read_all_vectors = True
    reader.read_all_tensors = True
    reader.read_all_fields = True
    err = reader.update()
    data = reader.get_output()
    
    
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
    if cellId == 5:
      ncorners = 3
    elif cellId == 3:
      ncorners = 2
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
