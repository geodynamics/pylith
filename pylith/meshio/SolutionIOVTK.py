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

## @file pyre/meshio/SolutionIOVTK.py
##
## @brief Python object for writing solution of finite-element problem
## to VTK file.
##
## Factory: solution_io

from SolutionIO import SolutionIO

# SolutionIOVTK class
class SolutionIOVTK(SolutionIO):
  """
  Python object for writing solution of finite-element problem to VTK file.

  Factory: solution_io
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(SolutionIO.Inventory):
    """
    Python object for managing SolutionIOVTK facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing SolutionIOVTK facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Name of mesh file
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="output.vtk")
    filename.meta['tip'] = "Name of VTK file."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solutioniovtk"):
    """
    Constructor.
    """
    SolutionIO.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    SolutionIO._configure(self)
    self.filename = self.inventory.filename
    return


  def _sync(self):
    """
    Force synchronization between Python and C++.
    """
    SolutionIO._sync(self)
    self.cppHandle.filename = self.filename
    return
  

  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.meshio.meshio as bindings
      self.cppHandle = bindings.SolutionIOVTK()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def solution_io():
  """
  Factory associated with SolutionIOVTK.
  """
  return SolutionIOVTK()


# End of file 
