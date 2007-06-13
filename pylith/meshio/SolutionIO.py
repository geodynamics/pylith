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
    ## @li \b output_freq Flag indicating whether to use 'time_step' or 'skip'
    ##   to set frequency of solution output.
    ## @li \b time_step Time step between solution output.
    ## @li \b skip Number of time steps to skip between solution output.
    ##
    ## \b Facilities
    ## @li \b coordsys Coordinate system for output.

    import pyre.inventory

    outputFreq = pyre.inventory.str("output_freq", default="skip",
             validator=pyre.inventory.choice(["skip", "time_step"]))
    outputFreq.meta['tip'] = "Flag indicating whether to use 'time_step' " \
                             "or 'skip' to set frequency of solution output."

    from pyre.units.time import s
    dt = pyre.inventory.dimensional("time_step", default=1.0*s)
    dt.meta['tip'] = "Time step between solution output."

    skip = pyre.inventory.int("skip", default=0,
                              validator=pyre.inventory.greaterEqual(0))
    skip.meta['tip'] = "Number of time steps to skip between solution output."

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
    self.mesh = None
    self.t = None
    self.istep = None
    return


  def open(self, mesh):
    """
    Open files for solution.
    """
    self._info.log("Opening files for output of solution.")
    self.mesh = mesh

    # Set flags
    self._sync()

    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for output is unknown."
    self.coordsys.initialize()

    assert(self.cppHandle != None)
    self.cppHandle.open(mesh.cppHandle)
    return


  def close(self):
    """
    Close files for solution.
    """
    self._info.log("Closing files for output of solution.")

    # Set flags
    self._sync()

    assert(self.cppHandle != None)
    self.cppHandle.close()
    return


  def writeTopology(self):
    """
    Write solution topology to file.
    """
    self._info.log("Writing solution topology.")

    assert(self.cppHandle != None)
    assert(self.mesh.cppHandle != None)
    assert(self.mesh.coordsys.cppHandle != None)
    self.cppHandle.writeTopology(self.mesh.cppHandle,
                                 self.mesh.coordsys.cppHandle)
    return


  def writeField(self, t, istep, field, name):
    """
    Write solution field at time t to file.
    """
    write = False
    if self.istep == None or not "value" in dir(self.t):
      write = True
    elif self.outputFreq == "skip":
      if istep > self.istep + self.skip:
        write = True
    elif t >= self.t + self.dt:
      write = True
    if write:
      self._info.log("Writing solution field '%s'." % name)
      assert(self.cppHandle != None)
      assert(self.mesh.cppHandle != None)
      self.cppHandle.writeField(t.value, field, name, self.mesh.cppHandle)
      self.istep = istep
      self.t = t
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.outputFreq = self.inventory.outputFreq
    self.dt = self.inventory.dt
    self.skip = self.inventory.skip
    self.coordsys = self.inventory.coordsys
    return


  def _sync(self):
    """
    Force synchronization between Python and C++.
    """
    self.cppHandle.coordsys = self.coordsys.cppHandle
    return


# End of file 
