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

## @file pylith/faults/FaultCohesive.py
##

## @brief Python abstract base class for a fault surface implemented
## with cohesive elements.
##
## Factory: fault

from Fault import Fault
from faults import FaultCohesive as ModuleFaultCohesive

# FaultCohesive class
class FaultCohesive(Fault, ModuleFaultCohesive):
  """
  Python abstract base class for a fault surface implemeted with
  cohesive elements.

  Inventory

  @class Inventory
  Python object for managing FaultCohesive facilities and properties.
  
  \b Properties
  @li \b use_fault_mesh If true, use fault mesh to define fault;
    otherwise, use group of vertices to define fault.
  
  \b Facilities
  @li \b fault_mesh_importer Importer for fault mesh.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  useMesh = pyre.inventory.bool("use_fault_mesh", default=False)
  useMesh.meta['tip'] = "If true, use fault mesh to define fault; " \
      "otherwise, use group of vertices to define fault."

  # Future, improved implementation
  #from pylith.meshio.MeshIOAscii imoport MeshIOAscii
  #faultMeshImporter = pyre.inventory.facility("fault_mesh_importer",
  #                                            factory=MeshIOLagrit,
  #                                            family="mesh_io")
  #faultMeshImporter.meta['tip'] = "Importer for fault mesh."

  # Current kludge
  meshFilename = pyre.inventory.str("fault_mesh_filename",
                                    default="fault.inp")
  meshFilename.meta['tip'] = "Filename for fault mesh UCD file."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesive"):
    """
    Constructor.
    """
    Fault.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Fault._configure(self)
    ModuleFaultCohesive.useFaultMesh(self, self.inventory.useMesh)
    #ModuleFaultCohesive.faultMeshImporter(self, 
    #                                      self.inventory.faultMeshImporter)

    # Hardwire collocated quadrature
    self.faultQuadrature.inventory.cell._configure()
    self.faultQuadrature._configure()
    self.faultQuadrature.cell.collocateQuad = True
    self.faultQuadrature.cell.order = self.faultQuadrature.cell.degree
    return

  
# End of file 
