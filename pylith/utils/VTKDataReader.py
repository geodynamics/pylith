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
    cells = data.get_cells().to_array()


    data = {'vertices': None,
            'cells': None,
            'vertex_fields': None,
            'cell_fields': None,
            }
    return data


# End of file 
