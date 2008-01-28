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

## @file pylith/faults/Fault.py
##

## @brief Python abstract base class for a fault surface.
##
## This implementation of a fault associates both physical
## properties and a quadrature scheme with the fault.
##
## Factory: fault

from pyre.components.Component import Component

# Validator for direction
def validateDir(value):
  """
  Validate direction.
  """
  msg = "Direction must be a 3 component vector (list)."
  if not isinstance(value, list):
    raise ValueError(msg)
  if 3 != len(value):
    raise ValueError(msg)
  try:
    nums = map(float, value)
  except:
    raise ValueError(msg)
  return value


# Fault class
class Fault(Component):
  """
  Python abstract base class for a fault surface.

  This implementation of a fault associates both physical
  properties and a quadrature scheme with the fault.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Fault facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Fault facilities and properties.
    ##
    ## \b Properties
    ## @li \b id Fault identifier
    ## @li \b name Name of fault
    ## @li \b up_dir Up-dip or up direction
    ##   (perpendicular to along-strike and not collinear with fault normal;
    ##   only applies to fault surfaces in a 3-D domain).
    ## @li \b normal_dir General preferred direction for fault normal
    ##   (used to pick which of two possible normal directions for
    ##   interface; only applies to fault surfaces in a 3-D domain).
    ## @li \b mat_db Spatial database for bulk material properties
    ##   (used in improving conditioning of Jacobian matrix).
    ##
    ## \b Facilities
    ## @li \b quadrature Quadrature object for numerical integration

    import pyre.inventory

    id = pyre.inventory.int("id", default=100)
    id.meta['tip'] = "Fault identifier (must be unique across all faults " \
                     "and materials)."

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Name of fault."

    upDir = pyre.inventory.list("up_dir", default=[0, 0, 1],
                                validator=validateDir)
    upDir.meta['tip'] = "Up-dip or up direction " \
                        "(perpendicular to along-strike and not collinear " \
                        "with fault normal; only applies to fault surface " \
                        "in a 3-D domain)."

    normalDir = pyre.inventory.list("normal_dir", default=[1, 0, 0],
                                validator=validateDir)
    normalDir.meta['tip'] = "General preferred direction for fault normal " \
                 "(used to pick which of two possible normal directions for " \
                 "interface; only applies to fault surfaces in a 3-D domain)."

    from pylith.feassemble.quadrature.Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    matDB = pyre.inventory.facility("mat_db", family="spatial_database",
                                   factory=SimpleDB, args=["bulk materials"])
    matDB.meta['tip'] = "Spatial database for bulk material properties " \
                        "(used in improving conditioning of Jacobian matrix)."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fault"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="fault")
    self.cppHandle = None
    self.mesh = None
    return


  def adjustTopology(self, mesh):
    """
    Adjust mesh topology for fault implementation.
    """
    self._createCppHandle()
    
    assert(None != self.cppHandle)
    self.cppHandle.id = self.id
    self.cppHandle.label = self.label
    self.cppHandle.adjustTopology(mesh.cppHandle)
    return


  def preinitialize(self, mesh):
    """
    Setup fault.
    """
    self._createCppHandle()
    self.cppHandle.id = self.id
    self.cppHandle.label = self.label

    self.mesh = mesh
    
    self.quadrature.preinitialize()    
    self.cppHandle.quadrature = self.quadrature.cppHandle
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    faultDim = self.mesh.dimension() - 1
    if faultDim != self.quadrature.cell.cellDim:
      raise ValueError, \
            "Quadrature is incompatible with fault surface.\n" \
            "Dimensions for quadrature: %d, dimensions of fault: %d" % \
            (self.quadrature.cell.cellDim, faultDim)

    # :TODO: Make sure mesh has group of vertices with label.
    return
  

  def initialize(self, totalTime, numTimeSteps):
    """
    Initialize fault.
    """
    self.quadrature.initialize()
    self.matDB.initialize()

    assert(None != self.cppHandle)
    self.cppHandle.initialize(self.mesh.cppHandle,
                              self.mesh.coordsys.cppHandle,
                              self.upDir, self.normalDir,
                              self.matDB.cppHandle)
    from pylith.topology.Mesh import Mesh
    self.faultMesh = Mesh()
    self.faultMesh.initialize(self.mesh.coordsys)
    self.cppHandle.faultMesh(self.faultMesh.cppHandle)
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return self.faultMesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    self.id = self.inventory.id
    self.label = self.inventory.label
    self.upDir = map(float, self.inventory.upDir)
    self.normalDir = map(float, self.inventory.normalDir)
    self.quadrature = self.inventory.quadrature
    self.matDB = self.inventory.matDB
    return

  
  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    raise NotImplementedError("Please implement _createCppHandle() in " \
                              "derived class.")
  
  
# End of file 
