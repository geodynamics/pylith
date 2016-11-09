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

## @file pylith/problems/Solution.py
##
## @brief Python solution field for problem.
##
## Factory: solution.

from pylith.utils.PetscComponent import PetscComponent

# Solution class =======================================================
class Solution(PetscComponent):
  """
  Python solution field for problem.

  Factory: solution.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing Solution facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Solution facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b subfields Subfields in solution.

    import pyre.inventory

    from SolnDisp import SolnDisp
    from .SolutionSubfield import subfieldFactory
    subfields = pyre.inventory.facilityArray("subfields", itemFactory=subfieldFactory, factory=SolnDisp)
    subfields.meta['tip'] = "Subfields in solution."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solution"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="solution")
    return


  def preinitialize(self, mesh, normalizer):
    """
    """
    from pylith.topology.Field import Field
    self.solution = Field(mesh)
    spaceDim = mesh.coordsys().spaceDim
    for subfield in self.subfields.components():
      subfield.initialize(normalizer, spaceDim)
      self.solution.subfieldAdd(subfield.label, subfield.ncomponents, subfield.vectorFieldType, subfield.scale)
    self.solution.subfieldsSetup()
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def initialize(self, constraints, integrators):
    """
    Setup solution field.
    """
    solution.setupSolnChart()

    for constraint in constraints.components():
      constraint.setConstraints(solution)
    solution.allocate()
    solution.zeroAll()
    for integrator in self.integrators.components():
      integrator.checkConstraints(solution)

    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.subfields = self.inventory.subfields
    return


# FACTORIES ////////////////////////////////////////////////////////////

def solution():
  """
  Factory associated with Solution.
  """
  return Solution()


# End of file
