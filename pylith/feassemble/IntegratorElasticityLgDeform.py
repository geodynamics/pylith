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

## @file pylith/feassemble/IntegratorElasticityLgDeform.py
##
## @brief Python object implementing sgeneral methods for time
## integration of the elasticity equation using finite-elements with
## large rogid body motions and small strains.
##
## Factory: integrator

from IntegratorElasticity import IntegratorElasticity

# IntegratorElasticityLgDeform class
class IntegratorElasticityLgDeform(IntegratorElasticity):
  """
  Python object implementing sgeneral methods for time integration of
  the elasticity equation using finite-elements with large rigid body
  motions and small strains.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="integratorelasticitylgdeform"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self)
    self.name = "Integrator ElasticityLgDeform"

    return


# End of file 
