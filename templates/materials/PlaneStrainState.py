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

## @brief Python object implementing 1-D isotropic linear elastic
## material for plane strain.
##
## Factory: material.

# ISA ElasticMaterial
from pylith.materials.ElasticMaterial import ElasticMaterial

# Import the SWIG module PlanseStrainState object and rename it
# ModulePlaneStrainState so that it doesn't clash with the local
# Python class with the same name.
from materialscontrib import PlaneStrainState as ModulePlaneStrainState

# PlaneStrainState class
class PlaneStrainState(ElasticMaterial, ModulePlaneStrainState):
  """
  Python object implementing 2-D isotropic linear elastic material for
  plane strain.

  Factory: material.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="planestrainstate"):
    """
    Constructor.
    """
    ElasticMaterial.__init__(self, name)
    # Set the fields that are available for output. These are the
    # stored physical properties, state variables, and the total
    # strain tensor and the stress tensor. For bulk elasticity
    # materials we can compute the stresses and strains in a general
    # fashion, so they need not be stored as they are in this example.
    #
    # There are no vertex fields because the constitutive model
    # operations on quantities evaluated at the quadrature points.
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': ["mu", "lambda", "density"],
            'data': ["total_strain", "stress"]}}
    self._loggingPrefix = "MaPlSn " # Prefix that appears in PETSc logging
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModulePlaneStrainState.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def material(): # The name of this function MUST be 'material'.
  """
  Factory associated with PlaneStrainState.
  """
  return PlaneStrainState() # Return our object


# End of file 
