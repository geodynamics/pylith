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

    includeDealloc = pyre.inventory.bool("include_dealloc", default=True)
    includeDealloc.meta['tip'] = "Subtract deallocated memory when reporting."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="perf_logger"):
    """
    Constructor.
    """
    Logger.__init__(self, name)
    self.megabyte = float(2**20)
    self.memory   = {}
    self.memory['Completion'] = 0
    return

  def logMesh(self, stage, mesh):
    """
    Read mesh parameters to determine memory from our model.
    """
    import pylith.perf.Mesh

    if not stage in self.memory: self.memory[stage] = {}
    meshModel = pylith.perf.Mesh.Mesh(mesh.dimension(), mesh.coneSize(), mesh.numVertices(), mesh.numCells())
    meshModel.tabulate(self.memory[stage])
    for group, nvertices in mesh.groupSizes():
      self.logVertexGroup('VertexGroups', group, nvertices, mesh.numVertices())
    return

  def logVertexGroup(self, stage, label, nvertices, nMeshVertices):
    """
    Read vertex group parameters to determine memory from our model.
    """
    import pylith.perf.VertexGroup

    if not stage in self.memory: self.memory[stage] = {}
    groupModel = pylith.perf.VertexGroup.VertexGroup(label, nvertices, nMeshVertices)
    groupModel.tabulate(self.memory[stage])
    return

  def logMaterial(self, stage, material):
    """
    Read material parameters to determine memory from our model.
    """
    import pylith.perf.Material

    if not stage in self.memory: self.memory[stage] = {}
    materialModel = pylith.perf.Material.Material(material.label(), material.ncells)
    materialModel.tabulate(self.memory[stage])
    self.logField(stage, material.propertiesField())
    self.logField(stage, material.stateVarsField())
    return

  def logQuadrature(self, stage, quadrature):
    ##self.logField(stage, quadrature.quadPtsPrecomp())
    ##self.logField(stage, quadrature.jacobianPrecomp())
    ##self.logField(stage, quadrature.jacobianDetPrecomp())
    ##self.logField(stage, quadrature.basisDerivPrecomp())
    # For debugging right now
    ##self.logField('Field', quadrature.quadPtsPrecomp())
    ##self.logField('Field', quadrature.jacobianPrecomp())
    ##self.logField('Field', quadrature.jacobianDetPrecomp())
    ##self.logField('Field', quadrature.basisDerivPrecomp())
    return

  def logField(self, stage, field):
    """
    Read field parameters to determine memory from our model.
    """
    import pylith.perf.Field

    if not stage in self.memory: self.memory[stage] = {}
    fieldModel = pylith.perf.Field.Field(field.label(), field.sectionSize(), field.chartSize())
    fieldModel.tabulate(self.memory[stage])
    return

  def logGlobalOrder(self, stage, label, field):
    """
    Read parameters to determine memory from our model.
    """
    import pylith.perf.GlobalOrder

    if not stage in self.memory: self.memory[stage] = {}
    orderModel = pylith.perf.GlobalOrder.GlobalOrder(label, field.chartSize())
    orderModel.tabulate(self.memory[stage])
    return

  def logJacobian(self, stage, label):
    """
    Read parameters to determine memory from our model.
    """
    import pylith.perf.Jacobian

    if not stage in self.memory: self.memory[stage] = {}
    jacModel = pylith.perf.Jacobian.Jacobian(label)
    jacModel.tabulate(self.memory[stage])
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

  def processMemDict(self, memDict, indent = 0, namePrefix = '', includeDealloc = True):
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
        out,mem,codeMem = self.processMemDict(m, indent, fullname, includeDealloc)
        output.extend(out)
        total     += mem
        codeTotal += codeMem
      else:
        mem = logger.getAllocationTotal(fullname)
        if includeDealloc:
          mem     -= logger.getDeallocationTotal(fullname)
        total     += m
        codeTotal += mem
        output.append(self.memLine('Model', name, m, indent))
        if logger.getAllocationTotal(fullname):
          output.append(self.memLine('Code',  name, mem, indent))
    if namePrefix:
      mem = logger.getAllocationTotal(namePrefix)
      if includeDealloc:
        mem -= logger.getDeallocationTotal(namePrefix)
    else:
      mem = logger.getAllocationTotal()
      if includeDealloc:
        mem -= logger.getDeallocationTotal()
    if mem == 0:
      mem = codeTotal
    output.append(self.memLine('Model', 'Total', total, indent))
    output.append(self.memLine('Code',  'Total', mem, indent))
    if mem:
      output.append('%sPercentage memory modeled: %.2f%%' % (self.prefix(indent), total*100.0/mem))
    else:
      output.append('%sMemory modeled:            %d (no measured memory)' % (self.prefix(indent), total))
    return output, total, codeTotal

  def show(self):
    """
    Print memory usage.
    """
    output = ["MEMORY USAGE"]
    output.extend(self.processMemDict(self.memory, includeDealloc = self.includeDealloc)[0])

    from pylith.utils.petsc import MemoryLogger
    logger    =  MemoryLogger.singleton()
    output.append(self.memLine('Code',  'Total Alloced', logger.getAllocationTotal('default'), 1))
    output.append(self.memLine('Code',  'Total Dealloced', logger.getDeallocationTotal('default'), 1))

    print '\n'.join(output)
    return

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Logger._configure(self)
    self.includeDealloc = self.inventory.includeDealloc
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
