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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

from pyre.components.Component import Component
from pyre.applications.Script import Script as Application

sizeInt = 4
sizeDouble = 8
import distutils.sysconfig

pointerSize = distutils.sysconfig.get_config_var('SIZEOF_VOID_P')
if pointerSize == 4:
  sizeArrow = 40 # 32 bit
elif pointerSize == 8:
  sizeArrow = 56 # 64 bit
else:
  raise RuntimeError('Could not determine the size of a pointer')

# ITEM FACTORIES ///////////////////////////////////////////////////////

def materialFactory(name):
  """
  Factory for material items.
  """
  from pyre.inventory import facility
  return facility(name, factory=Material)


def vertexGroupFactory(name):
  """
  Factory for vertex group items.
  """
  from pyre.inventory import facility
  return facility(name, factory=VertexGroup)


# ----------------------------------------------------------------------
# Mesh class
import pylith.perf.Mesh
class Mesh(Component,pylith.perf.Mesh.Mesh):
  """
  Mesh object for holding mesh size information.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Mesh facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Mesh facilities and properties.
    ##
    ## \b Properties
    ## @li \b cell_type Type of cell
    ## @li \b ncells Number of cells
    ## @li \b nvertices Number of vertices
    ## @li \b ngroups Number of vertex groups
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    cellType = pyre.inventory.str('cell_type', default="hex8",
         validator=pyre.inventory.choice(['tri3', 'quad4', 'tet4', 'hex8']))
    cellType.meta['tip'] = "Type of cell in mesh."

    ncells = pyre.inventory.int('ncells', default=0)
    ncells.meta['tip'] = "Number of cells in finite-element mesh."

    nvertices = pyre.inventory.int('nvertices', default=0)
    nvertices.meta['tip'] = "Number of vertices in finite-element mesh."

    ngroups = pyre.inventory.int('ngroups', default=0)
    ngroups.meta['tip'] = "Number of vertex groups in finite-element mesh."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="mesh"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    dim, coneSize = pylith.perf.Mesh.Mesh.cellTypeInfo(self.inventory.cellType)
    pylith.perf.Mesh.Mesh.__init__(self, dim, coneSize, self.inventory.nvertices, self.inventory.ncells)
    return


# ----------------------------------------------------------------------
# VertexGroup class
import pylith.perf.VertexGroup
class VertexGroup(Component,pylith.perf.VertexGroup.VertexGroup):
  """
  Mesh object for holding vertex group size information.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing VertexGroup facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing VertexGroup facilities and properties.
    ##
    ## \b Properties
    ## @li \b label Name of group.
    ## @li \b size Number of vertices in group.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    label = pyre.inventory.str('label', default="Vertex Group")
    label.meta['tip'] = "Label for vertex group."

    size = pyre.inventory.int('size', default=0)
    size.meta['tip'] = "Number of vertices in group."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="vertexgroup"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    pylith.perf.VertexGroup.VertexGroup.__init__(self, self.inventory.label, self.inventory.size)
    return


# ----------------------------------------------------------------------
# Material class
import pylith.perf.Material
class Material(Component,pylith.perf.Material.Material):
  """
  Mesh object for holding material size information.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Material facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Mesh facilities and properties.
    ##
    ## \b Properties
    ## @li \b label Name of material.
    ## @li \b size Number of cells in material.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    label = pyre.inventory.str('label', default="Material")
    label.meta['tip'] = "Label for material."

    size = pyre.inventory.int('size', default=0)
    size.meta['tip'] = "Number of cells in material."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="material"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    pylith.perf.Material.Material.__init__(self, self.inventory.label, self.inventory.size)
    return


# ----------------------------------------------------------------------
# MemoryUsageApp class
class MemoryUsageApp(Application):
  """
  Python MemoryUsage application.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Application.Inventory):
    """
    Python object for managing MemoryUsageApp facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MemoryUsageApp facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li mesh Finite-element mesh information.
    ## @li groups Vertex groups information.
    ## @li materials Materials information.

    import pyre.inventory
    from pylith.utils.EmptyBin import EmptyBin

    mesh = pyre.inventory.facility("mesh", factory=Mesh)
    mesh.meta['tip'] = "Finite-element mesh information."

    groups = pyre.inventory.facilityArray("vertex_groups",
                                          itemFactory=vertexGroupFactory,
                                          factory=EmptyBin)
    groups.meta['tip'] = "Vertex groups information."

    materials = pyre.inventory.facilityArray("materials",
                                          itemFactory=materialFactory,
                                             factory=EmptyBin)
    materials.meta['tip'] = "Materials information."

    
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="memoryusage"):
    """
    Constructor.
    """
    Application.__init__(self, name)
    self.memory = {'mesh': 0,
                   'stratification': 0,
                   'groups': 0}
    self.mesh = None
    self.groups = None
    self.materials = None
    return


  def main(self, *args, **kwds):
    """
    Run the application.
    """
    self.mesh.initialize()
    self._tabulate()
    self.show()
    return
  

  def show(self):
    """
    Print memory usage.
    """
    megabyte = float(2**20)

    print "MEMORY USAGE"
    print "  Finite-element mesh"
    memory = self.memory['mesh']
    print "    Mesh:           %d bytes (%.3f MB)" % \
          (memory['mesh'],
           memory['mesh'] / megabyte)
    print "    Stratification: %d bytes (%.3f MB)" % \
          (memory['stratification'],
           memory['stratification'] / megabyte)
    print "    Coordinates:    %d bytes (%.3f MB)" % \
          (memory['coordinates'],
           memory['coordinates'] / megabyte)
    print "    Materials:      %d bytes (%.3f MB)" % \
          (memory['materials'],
           memory['materials'] / megabyte)

    print "    Groups"
    for (label, nbytes) in self.memory['groups'].items():
      print "      %s: %d bytes (%.3f MB)" % \
            (label, nbytes, nbytes / megabyte)

    print "  Materials"
    for (label, nbytes) in self.memory['materials'].items():
      print "    %s: %d bytes (%.3f MB)" % \
            (label, nbytes, nbytes / megabyte)
    print "  TOTAL:            %d bytes (%.3f MB)" % \
          (self.memory['total'], self.memory['total'] / megabyte)

    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.mesh = self.inventory.mesh
    self.groups = self.inventory.groups
    self.materials = self.inventory.materials
    return


  def _tabulate(self):
    """
    Tabulate expected memory usage.
    """
    total = 0
    memory = {}

    # mesh
    info = self.mesh.tabulate()
    memory['mesh'] = info
    total += sum(info.values())

    # groups
    memory['groups'] = {}
    for group in self.groups.components():
      nbytes = group.tabulate()
      total += nbytes
      memory['groups'][group.label] = nbytes

    # materials
    memory['materials'] = {}
    for material in self.materials.components():
      nbytes = material.tabulate()
      total += nbytes
      memory['materials'][material.label] = nbytes

    memory['total'] = total
    self.memory = memory
    return
  

# ======================================================================
__requires__ = ""

if __name__ == "__main__":

  app = MemoryUsageApp()
  app.run()


# End of file 
