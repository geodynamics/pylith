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

  class Inventory(SolutionIO.Inventory):
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


  def write(self, t, fields, mesh):
    """
    Write solution at time t to file.
    """
    self._info.log("Writing solution information.")

    # Set flags
    self._sync()

    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for mesh is unknown."
    self.coordsys.initialize()

    #self.cppHandle.write(t, fields, mesh.cppHandle)
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
    #self.cppHandle.coordsys = self.coordsys
    return


# End of file 
