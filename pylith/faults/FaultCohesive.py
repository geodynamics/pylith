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

## @file pylith/faults/FaultCohesive.py
##

## @brief Python abstract base class for a fault surface implemented
## with cohesive elements.
##
## Factory: fault

from Fault import Fault

# FaultCohesive class
class FaultCohesive(Fault):
  """
  Python abstract base class for a fault surface implemeted with
  cohesive elements.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Fault.Inventory):
    """
    Python object for managing FaultCohesive facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FaultCohesive facilities and properties.
    ##
    ## \b Properties
    ## @li \b use_fault_mesh If true, use fault mesh to define fault;
    ##   otherwise, use group of vertices to define fault.
    ##
    ## \b Facilities
    ## @li \b fault_mesh_importer Importer for fault mesh.

    import pyre.inventory

    useFaultMesh = pyre.inventory.bool("use_fault_mesh", default=False)
    useFaultMesh.meta['tip'] = "If true, use fault mesh to define fault; " \
        "otherwise, use group of vertices to define fault."

    # Future, improved implementation
    #from pylith.meshio.MeshIOAscii imoport MeshIOAscii
    #faultMeshImporter = pyre.inventory.facility("fault_mesh_importer",
    #                                            factory=MeshIOLagrit,
    #                                            family="mesh_io")
    #faultMeshImporter.meta['tip'] = "Importer for fault mesh."

    # Current kludge
    faultMeshFilename = pyre.inventory.str("fault_mesh_filename",
                                           default="fault.inp")
    faultMeshFilename.meta['tip'] = "Filename for fault mesh UCD file."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesive"):
    """
    Constructor.
    """
    Fault.__init__(self, name)
    return


  def adjustTopology(self, mesh):
    """
    Adjust mesh topology for fault implementation.
    """
    self._createCppHandle()
    assert(None != self.cppHandle)
    self.cppHandle.useFaultMesh = self.useFaultMesh
    #self.cppHandle.faultMeshImporter = self.faultMeshImporter.cppHandle
    self.cppHandle.faultMeshFilename = self.faultMeshFilename # TEMPORARY

    Fault.adjustTopology(self, mesh)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Fault._configure(self)
    self.useFaultMesh = self.inventory.useFaultMesh
    #self.faultMeshImporter = self.inventory.faultMeshImporter
    self.faultMeshFilename = self.inventory.faultMeshFilename # TEMPORARY
    return

  
# End of file 
