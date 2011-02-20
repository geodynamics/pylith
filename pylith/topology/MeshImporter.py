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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/topology/MeshImporter.py
##
## @brief Python implementation of importing a mesh.
##
## Factory: mesh_generator.

from MeshGenerator import MeshGenerator

# MeshImporter class
class MeshImporter(MeshGenerator):
  """
  Python implementation of importing a mesh.

  Factory: mesh_generator.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MeshGenerator.Inventory):
    """
    Python object for managing MeshImporter facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MeshImporter facilities and properties.
    ##
    ## \b Properties
    ## @li reorder_mesh Reorder mesh using reverse Cuthill-McKee if true.
    ##
    ## \b Facilities
    ## @li \b reader Mesh reader.
    ## @li \b partitioner Mesh partitioner.
    ## @li \b refiner Mesh refiner.

    import pyre.inventory

    reorderMesh = pyre.inventory.bool("reorder_mesh", default=False)
    reorderMesh.meta['tip'] = "Reorder mesh using reverse Cuthill-McKee."

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    reader = pyre.inventory.facility("reader", family="mesh_io",
                                       factory=MeshIOAscii)
    reader.meta['tip'] = "Mesh reader."

    from Distributor import Distributor
    distributor = pyre.inventory.facility("distributor",
                                          family="mesh_distributor",
                                          factory=Distributor)
    distributor.meta['tip'] = "Mesh distributor."
 
    from MeshRefiner import MeshRefiner
    refiner = pyre.inventory.facility("refiner",
                                      family="mesh_refiner",
                                      factory=MeshRefiner)
    refiner.meta['tip'] = "Mesh refiner."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshimporter"):
    """
    Constructor.
    """
    MeshGenerator.__init__(self, name)
    self._loggingPrefix = "MeIm "
    return


  def create(self, normalizer, faults=None):
    """
    Hook for creating mesh.
    """
    from pylith.utils.profiling import resourceUsageString

    self._setupLogging()
    logEvent = "%screate" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)    

    # Read mesh
    mesh = self.reader.read(self.debug, self.interpolate)
    if self.debug:
      mesh.view("Finite-element mesh.")

    # Reorder mesh
    if self.reorderMesh:
      logEvent2 = "%sreorder" % self._loggingPrefix
      self._eventLogger.eventBegin(logEvent2)
      self._debug.log(resourceUsageString())
      self._info.log("Reordering cells and vertices.")
      from pylith.topology.ReverseCuthillMcKee import ReverseCuthillMcKee
      ordering = ReverseCuthillMcKee()
      ordering.reorder(mesh)
      self._eventLogger.eventEnd(logEvent2)

    # Adjust topology
    self._debug.log(resourceUsageString())
    self._info.log("Adjusting topology.")
    self._adjustTopology(mesh, faults)

    # Distribute mesh
    if mesh.getComm().size > 1:
      self._info.log("Distributing mesh.")
      mesh = self.distributor.distribute(mesh, normalizer)
      if self.debug:
        mesh.view("Distributed mesh.")

    # Refine mesh (if necessary)
    mesh = self.refiner.refine(mesh)

    # Nondimensionalize mesh (coordinates of vertices).
    mesh.nondimensionalize(normalizer)

    self._eventLogger.eventEnd(logEvent)    
    return mesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based on inventory.
    """
    MeshGenerator._configure(self)
    self.reader = self.inventory.reader
    self.distributor = self.inventory.distributor
    self.refiner = self.inventory.refiner
    self.reorderMesh = self.inventory.reorderMesh
    return
  

  def _setupLogging(self):
    """
    Setup event logging.
    """
    MeshGenerator._setupLogging(self)
    self._eventLogger.registerEvent("%sreorder" % self._loggingPrefix)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_generator():
  """
  Factory associated with MeshImporter.
  """
  return MeshImporter()


# End of file 
