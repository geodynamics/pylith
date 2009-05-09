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
    self.megabyte = float(2**20)
    self.memory   = {}
    return

  def logMesh(self, stage, mesh):
    """
    Read mesh parameters to determine memory from our model.
    """
    import pylith.perf.Mesh

    if not 'Mesh' in self.memory: self.memory['Mesh'] = {}
    meshModel = pylith.perf.Mesh.Mesh(mesh.dimension(), mesh.coneSize(), mesh.numVertices(), mesh.numCells())
    meshModel.tabulate(self.memory['Mesh'])
    for group, nvertices in mesh.groupSizes():
      self.logVertexGroup(stage, group, nvertices, mesh.numVertices())
    if self.verbose: self.show()
    return

  def logVertexGroup(self, stage, label, nvertices, nMeshVertices):
    """
    Read vertex group parameters to determine memory from our model.
    """
    import pylith.perf.VertexGroup

    if not 'VertexGroups' in self.memory: self.memory['VertexGroups'] = {}
    group = pylith.perf.VertexGroup.VertexGroup(label, nvertices, nMeshVertices)
    group.tabulate(self.memory['VertexGroups'])
    if self.verbose: self.show()
    return

  def logMaterial(self, stage, material):
    """
    Read material parameters to determine memory from our model.
    """
    import pylith.perf.Material

    if not 'Materials' in self.memory: self.memory['Materials'] = {}
    material = pylith.perf.Material.Material(material.label(), material.ncells)
    material.tabulate(self.memory['Materials'])
    if self.verbose: self.show()
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

  def join(self, logger):
    """
    Incorporate information from another logger.
    """
    self.mergeMemDict(self.memory, logger.memory)
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

  def show(self):
    """
    Print memory usage.
    """
    output = ["MEMORY USAGE"]
    output.extend(self.processMemDict(self.memory)[0])
    print '\n'.join(output)
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
