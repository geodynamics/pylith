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

## @file pyre/meshio/SolutionIO.py
##
## @brief Python abstract base class for I/O of the finite-element
## solution.
##
## Factory: solution_io

from pyre.components.Component import Component

# SolutionIO class
class SolutionIO(Component):
  """
  Python abstract base class for finite-element mesh I/O.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing SolutionIOVTK facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing SolutionIOVTK facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b coordsys Coordinate system for output.

    import pyre.inventory

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system for output."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solutionio"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="solution_io")
    self.cppHandle = None
    self.coordsys = None
    return


  def open(self, mesh):
    """
    Open files for solution.
    """
    self._info.log("Opening files for output of solution.")

    # Set flags
    self._sync()

    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for output is unknown."
    self.coordsys.initialize()

    assert(cppHandle != None)
    self.cppHandle.open(mesh.cppHandle)
    return


  def writeTopology(self, mesh):
    """
    Write solution topology to file.
    """
    self._info.log("Writing solution topology.")

    assert(cppHandle != None)
    self.cppHandle.write(mesh.cppHandle, mesh.coordsys.cppHandle)
    return


  def writeField(self, t, field, name, mesh):
    """
    Write solution field at time t to file.
    """
    self._info.log("Writing solution field '%s'." % name)

    self.cppHandle.write(t, field, name, mesh.cppHandle)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.coordsys = self.inventory.coordsys
    return


  def _sync(self):
    """
    Force synchronization between Python and C++.
    """
    self.cppHandle.coordsys = self.coordsys
    return


# End of file 
