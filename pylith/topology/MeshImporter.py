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
    ## @li None
    ##
    ## \b Facilities
    ## @li \b reader Mesh reader.
    ## @li \b partitioner Mesh partitioner.
    ## @li \b refiner Mesh refiner.

    import pyre.inventory

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
    self._logger.eventBegin(logEvent)    

    mesh = self.reader.read(normalizer, self.debug, self.interpolate)
    if self.debug:
      mesh.view("Finite-element mesh.")
    self._debug.log(resourceUsageString())
    self._info.log("Adjusting topology.")
    self._adjustTopology(mesh, faults)

    import mpi
    if mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD) > 1:
      self._info.log("Distributing mesh.")
      mesh = self.distributor.distribute(mesh, normalizer)
      if self.debug:
        mesh.view("Distributed mesh.")

    # refine mesh (if necessary)
    mesh = self.refiner.refine(mesh)

    self._logger.eventEnd(logEvent)    
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
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_generator():
  """
  Factory associated with MeshImporter.
  """
  return MeshImporter()


# End of file 
