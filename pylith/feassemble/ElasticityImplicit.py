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

## @file pylith/feassemble/ElasticityImplicit.py
##
## @brief Python object for implicit time integration of dynamic
## elasticity equation using finite-elements.
##
## Factory: integrator

from IntegratorElasticity import IntegratorElasticity

# ElasticityImplicit class
class ElasticityImplicit(IntegratorElasticity):
  """
  Python object for implicit time integration of dynamic elasticity
  equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticityimplicit"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self, name)

    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.ElasticityImplicit()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ElasticityImplicit.
  """
  return ElasticityImplicit()


# End of file 
