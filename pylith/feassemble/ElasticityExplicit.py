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

## @file pylith/feassemble/ElasticityExplicit.py
##
## @brief Python object for explicit time integration of dynamic
## elasticity equation using finite-elements.
##
## Factory: integrator

from IntegratorElasticity import IntegratorElasticity

# ElasticityExplicit class
class ElasticityExplicit(IntegratorElasticity):
  """
  Python object for explicit time integration of dynamic elasticity
  equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticityexplicit"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self, name)

    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.ElasticityExplicit()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ElasticityExplicit.
  """
  return ElasticityExplicit()


# End of file 
