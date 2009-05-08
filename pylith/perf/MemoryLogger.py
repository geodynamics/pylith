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

## @file pylith/perf/MemoryLogger.py
##
## @brief Python base class for performance and memory logging.
##
## Factory: perf_logger.

from Logger import Logger

# MemoryLogger class
class MemoryLogger(Logger):
  """
  Python base class for performance and memory logging.

  Factory: perf_logger.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Logger.Inventory):
    """
    Python object for managing Problem facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Problem facilities and properties.
    ##
    ## \b Properties
    ## @li \b dummy Nothing.

    import pyre.inventory

    dummy = pyre.inventory.bool("dummy", default=True)
    dummy.meta['tip'] = "Nothing."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="perf_logger"):
    """
    Constructor.
    """
    Logger.__init__(self, name)
    self.megabyte     = float(2**20)
    self.materials    = {}
    self.vertexGroups = {}
    return

  def logMesh(self, stage, mesh):
    """
    Read mesh parameters to determine memory from our model.
    """
    import pylith.perf.Mesh
    self.mesh = pylith.perf.Mesh.Mesh(mesh.dimension(), mesh.coneSize(), mesh.numVertices(), mesh.numCells())
    for group, nvertices in mesh.groupSizes():
      self.logVertexGroup(stage, group, nvertices)
    if self.verbose:
      self.tabulate()
      self.show()
    return

  def logVertexGroup(self, stage, label, nvertices):
    """
    Read vertex group parameters to determine memory from our model.
    """
    import pylith.perf.VertexGroup
    self.vertexGroups[label] = pylith.perf.VertexGroup.VertexGroup(label, nvertices)
    if self.verbose:
      self.tabulate()
      self.show()
    return

  def logMaterial(self, stage, material):
    """
    Read material parameters to determine memory from our model.
    """
    import pylith.perf.Material
    self.materials[material.id()] = pylith.perf.Material.Material(material.label(), material.ncells)
    if self.verbose:
      self.tabulate()
      self.show()
    return

  def tabulate(self):
    """
    Tabulate expected memory usage.
    """
    total = 0
    memory = {}
    # mesh
    if hasattr(self, 'mesh'):
      info = self.mesh.tabulate()
      memory['mesh'] = info
      total += sum(info.values())
    # groups
    if hasattr(self, 'vertexGroups'):
      memory['groups'] = {}
      for label,group in self.vertexGroups.iteritems():
        nbytes = group.tabulate()
        total += nbytes
        memory['groups'][label] = nbytes
    # materials
    if hasattr(self, 'materials'):
      memory['materials'] = {}
      for id, material in self.materials.iteritems():
        nbytes = material.tabulate()
        total += nbytes
        memory['materials'][material.label] = nbytes
    # Total
    memory['total'] = total
    self.memory = memory
    return

  def join(self, logger):
    """
    Incorporate information from another logger.
    """
    # mesh
    memory = logger.memory
    if 'mesh' in memory:
      if not 'mesh' in self.memory:
        self.memory['mesh'] = memory['mesh']
      else:
        for key in self.memory['mesh']:
          self.memory['mesh'][key] += memory['mesh'][key]
    # groups
    if 'groups' in memory:
      if not 'groups' in self.memory:
        self.memory['groups'] = memory['groups']
      else:
        for key in memory['groups']:
          if not key in self.memory['groups']:
            self.memory['groups'][key]  = memory['groups'][key]
          else:
            self.memory['groups'][key] += memory['groups'][key]
    # materials
    if 'materials' in memory:
      if not 'materials' in self.memory:
        self.memory['materials'] = memory['materials']
      else:
        for key in memory['materials']:
          if not key in self.memory['materials']:
            self.memory['materials'][key]  = memory['materials'][key]
          else:
            self.memory['materials'][key] += memory['materials'][key]
    # Total
    self.memory['total'] += memory['total']
    return

  def show(self):
    """
    Print memory usage.
    """
    megabyte = self.megabyte
    print "MEMORY USAGE"
    if 'mesh' in self.memory:
      print "  Finite-element mesh"
      memory = self.memory['mesh']
      print "    Mesh:           %d bytes (%.3f MB)" % (memory['mesh'], memory['mesh'] / megabyte)
      print "    Stratification: %d bytes (%.3f MB)" % (memory['stratification'], memory['stratification'] / megabyte)
      print "    Coordinates:    %d bytes (%.3f MB)" % (memory['coordinates'], memory['coordinates'] / megabyte)
      print "    Materials:      %d bytes (%.3f MB)" % (memory['materials'], memory['materials'] / megabyte)
    if 'groups' in self.memory:
      print "    Groups"
      for (label, nbytes) in self.memory['groups'].items():
        print "      %s: %d bytes (%.3f MB)" % (label, nbytes, nbytes / megabyte)
    if 'materials' in self.memory:
      print "  Materials"
      for (label, nbytes) in self.memory['materials'].items():
        print "    %s: %d bytes (%.3f MB)" % (label, nbytes, nbytes / megabyte)
    print "  TOTAL:           %d bytes (%.3f MB)" % (self.memory['total'], self.memory['total'] / megabyte)
    return

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Logger._configure(self)
    self.dummy = self.inventory.dummy
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def perf_logger():
  """
  Factory associated with Logger.
  """
  return MemoryLogger()


# End of file 
