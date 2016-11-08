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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/materials/IsotropicLinearElasticityPlaneStrain.py
##
## @brief Python material for isotropic, linearly elastic, plane
## strain material.
##
## Factory: material

from .MaterialNew import MaterialNew
from .materials import IsotropicLinearElasticityPlaneStrain as ModuleMaterial

# IsotropicLinearElasticityPlaneStrain class
class IsotropicLinearElasticityPlaneStrain(MaterialNew, ModuleMaterial):
  """
  Python material for isotropic, linearly elastic plane strain.

  Factory: material
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MaterialNew.Inventory):
    """
    Python object for managing IsotropicLinearElasticityPlaneStrain facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Material facilities and properties.
    ##
    ## \b Properties
    ## @li \b id Material identifier (from mesh generator)
    ## @li \b label Descriptive label for material.
    ##
    ## \b Facilities
    ## @li \b db_properties Database of material property parameters
    ## @li \b quadrature Quadrature object for numerical integration
    ## @li \b db_initial_state Database for initial state.

    import pyre.inventory

    useInertia = pyre.inventory.bool("use_inertia", default=False)
    useInertia.meta['tip'] = "Include inertial term in elasticity equation."

    useBodyForce = pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in elasticity equation."

    useReferenceState = pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    from .AuxFieldsIsotropicLinearElasticity import AuxFieldsIsotropicLinearElasticity
    from pylith.topology.AuxSubfield import subfieldFactory
    auxFields = pyre.inventory.facilityArray("auxiliary_fields", itemFactory=subfieldFactory, factory=AuxFieldsIsotropicLinearElasticity)
    auxFields.meta['tip'] = "Discretization of physical properties and state variables."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="isotropiclinearelasticityplanestrain"):
    """
    Constructor.
    """
    MaterialNew.__init__(self, name)
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    MaterialNew.preinitialize(self, mesh)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    try:
      MaterialNew._configure(self)
      ModuleMaterial.useInertia(self, self.inventory.useInertia)
      ModuleMaterial.useBodyForce(self, self.inventory.useBodyForce)
      ModuleMaterial.useReferenceState(self, self.inventory.useReferenceState)

    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring material "
                       "(%s):\n%s" % (aliases, err.message))
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleMaterial.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def material():
  """
  Factory associated with IsotropicLinearElasticityPlaneStrain.
  """
  return IsotropicLinearElasticityPlaneStrain()


# End of file
