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
    reader.file_name = "axialplanestrain-statevars-elastic_info.vtk"
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
    id = cellTypes[0]
    if numpy.sum(cellTypes-id) != 0:
      raise ValueError("Expecting cells to all be the same type.")
    if id == 5:
      ncorners = 3
    else:
      raise ValueError("Unknown VTK cell type '%d'." % id)
    cells = cells.reshape( (ncells, 1+ncorners) )[:,1:1+ncorners]
    return cells


  def _getVertexFields(self, data):
    return


  def _getCellFields(self, data):
    return


# End of file 
