# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

# @brief Python object implementing 1-D isotropic linear elastic
# material for plane strain.
##
# Factory: material.

# ISA ElasticMaterial
from pylith.materials.ElasticMaterial import ElasticMaterial

# Import the SWIG module PlanseStrainState object and rename it
# ModulePlaneStrainState so that it doesn't clash with the local
# Python class of the same name.
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
        #
        # Do not change the name of this variable. The output manager will
        # request this variable by name.
        self.availableFields = \
            {'vertex':
             {'info': [],
              'data': []},
             'cell':
             {'info': ["mu", "lambda", "density"],
              'data': ["total_strain", "stress"]}}
        self._loggingPrefix = "MaPlSn "  # Prefix that appears in PETSc logging
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object. This
        function is called automatically by the generic Python Material
        object. It must have this name and self as the only argument.
        """
        ModulePlaneStrainState.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

# This is the function that is called when you invoke
# material_one = pylith.materials.contrib.PlaneStrainState
# The name of this function MUST be 'material'.
def material():
    """
    Factory associated with PlaneStrainState.
    """
    return PlaneStrainState()  # Return our object


# End of file
