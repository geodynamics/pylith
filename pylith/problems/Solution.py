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

# OneField class
class OneField(PetscComponent):
  """
  Python subfields container with one subfield.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing on OneField facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Homogeneous facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b subfield Subfield in solution.

    import pyre.inventory

    from SolutionSubfield import SolutionSubfield
    subfield = pyre.inventory.facility("subfield", family="subfield", factory=Subfield)
    subfield.meta['tip'] = "Subfield in solution."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="onefield"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="onefield")
    return


# ITEM FACTORIES ///////////////////////////////////////////////////////

def subfieldFactory(name):
  """
  Factory for subfield items.
  """
  from pyre.inventory import facility
  from pylith.problems.SolutionSubfield import SolutionSubfield
  return facility(name, family="subfield", factory=SolutionSubfield)


# Solution class
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

    subfields = pyre.inventory.facilityArray("subfields", itemFactory=subfieldFactory, factory=OneField)
    subfields.meta['tip'] = "Subfields in solution."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solution"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="solution")
    return


  def preinitialize(self, mesh):
    """
    """
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
    solution = Field()
    for subfield in self.subfields:
      solution.subfieldAdd(subfield.label, subfield.ncomponents, subfield.vectorFieldType, subfield.scale)
    solution.subfieldsSetup()
    solution.setupSolnChart()

    for constraint in constraints.components():
      constraint.setConstraintSizes(solution)
    solution.allocate()
    solution.zeroAll()
    for constraint in self.constraints.components():
      constraint.setConstraints(solution)
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
