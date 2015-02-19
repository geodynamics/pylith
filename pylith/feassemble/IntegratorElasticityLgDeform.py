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
