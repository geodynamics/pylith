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
    self.memoryNew    = {}
    self.materials    = {}
    self.vertexGroups = {}
    return

  def logMesh(self, stage, mesh):
    """
    Read mesh parameters to determine memory from our model.
    """
    import pylith.perf.Mesh
    self.mesh = pylith.perf.Mesh.Mesh(mesh.dimension(), mesh.coneSize(), mesh.numVertices(), mesh.numCells())
    if not 'Mesh' in self.memoryNew: self.memoryNew['Mesh'] = {}
    self.mesh.tabulateNew(self.memoryNew['Mesh'])
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
    if not 'VertexGroups' in self.memoryNew: self.memoryNew['VertexGroups'] = {}
    self.vertexGroups[label].tabulateNew(self.memoryNew['VertexGroups'])
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
    if not 'Materials' in self.memoryNew: self.memoryNew['Materials'] = {}
    self.materials[material.id()].tabulateNew(self.memoryNew['Materials'])
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

  def mergeMemDict(self, memDictTarget, memDictSource):
    for key in memDictSource:
      if not key in memDictTarget:
        memDictTarget[key] = memDictSource[key]
      elif isinstance(memDictSource[key], dict):
        if not isinstance(memDictTarget[key], dict):
          raise RuntimeError('Type mismatch in memory dict for key '+str(key))
        self.mergeMemDict(memDictTarget[key], memDictSource[key])
      else:
        memDictTarget[key] += memDictSource[key]
    return

  def joinNew(self, logger):
    """
    Incorporate information from another logger.
    """
    self.mergeMemDict(self.memoryNew, logger.memoryNew)
    return

  def join(self, logger):
    """
    Incorporate information from another logger.
    """
    self.joinNew(logger)
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

  def prefix(self, indent):
    prefix = ''
    for i in range(indent):
      prefix += '  '
    return prefix

  def memLine(self, source, name, mem, indent = 0):
    return '%s%-30s %8d bytes (%.3f MB)' % (self.prefix(indent), name+' ('+source+'):', mem, mem / self.megabyte)

  def processMemDict(self, memDict, indent = 0, namePrefix = ''):
    from pylith.utils.petsc import MemoryLogger
    logger    =  MemoryLogger.singleton()
    output    = []
    total     = 0
    codeTotal = 0
    indent   += 1
    for name,m in memDict.iteritems():
      fullname = namePrefix+name
      if isinstance(m, dict):
        output.append(self.prefix(indent)+name)
        out,mem,codeMem = self.processMemDict(m, indent, fullname)
        output.extend(out)
        total     += mem
        codeTotal += codeMem
      else:
        mem = logger.getAllocationTotal(fullname) - logger.getDeallocationTotal(fullname)
        total     += m
        codeTotal += mem
        output.append(self.memLine('Model', name, m, indent))
        if logger.getAllocationTotal(fullname):
          output.append(self.memLine('Code',  name, mem, indent))
    if namePrefix:
      mem = logger.getAllocationTotal(namePrefix) - logger.getDeallocationTotal(namePrefix)
    else:
      mem = logger.getAllocationTotal() - logger.getDeallocationTotal()
    if mem == 0:
      mem = codeTotal
    output.append(self.memLine('Model', 'Total', total, indent))
    output.append(self.memLine('Code',  'Total', mem, indent))
    output.append('%sPercentage memory modeled: %.2f%%' % (self.prefix(indent), total*100.0/mem))
    return output, total, codeTotal

  def showNew(self):
    """
    Print memory usage.
    """
    output = ["MEMORY USAGE"]
    output.extend(self.processMemDict(self.memoryNew)[0])
    print '\n'.join(output)
    return

  def show(self):
    """
    Print memory usage.
    """
    from pylith.utils.petsc import MemoryLogger
    logger   = MemoryLogger.singleton()
    megabyte = self.megabyte
    print "MEMORY USAGE"
    if 'mesh' in self.memory:
      print "  Finite-element mesh"
      memory = self.memory['mesh']
      print "    Mesh:           %d bytes (%.3f MB)" % (memory['mesh'], memory['mesh'] / megabyte)
      print "    Stratification: %d bytes (%.3f MB)" % (memory['stratification'], memory['stratification'] / megabyte)
      print "    Coordinates:    %d bytes (%.3f MB)" % (memory['coordinates'], memory['coordinates'] / megabyte)
      print "    Materials:      %d bytes (%.3f MB)" % (memory['materials'], memory['materials'] / megabyte)
      stage  = "MeshCreation"
      mem    = logger.getAllocationTotal(stage) - logger.getDeallocationTotal(stage)
      print "    Mesh (Code):    %d bytes (%.3f MB)" % (mem, mem / megabyte)
    if 'groups' in self.memory:
      print "    Groups"
      for (label, nbytes) in self.memory['groups'].items():
        print "      %s: %d bytes (%.3f MB)" % (label, nbytes, nbytes / megabyte)
    if 'materials' in self.memory:
      print "  Materials"
      for (label, nbytes) in self.memory['materials'].items():
        print "    %s: %d bytes (%.3f MB)" % (label, nbytes, nbytes / megabyte)
    print "  TOTAL:           %d bytes (%.3f MB)" % (self.memory['total'], self.memory['total'] / megabyte)
    mem = logger.getAllocationTotal() - logger.getDeallocationTotal()
    print "  TOTAL (Code):    %d bytes (%.3f MB)" % (mem, mem / megabyte)
    print "  Percentage memory modeled: %.2f%%" % (self.memory['total']*100.0/mem)
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
