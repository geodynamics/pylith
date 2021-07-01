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

# @file pylith/friction/ViscousFriction.py
##
# @brief Python object implementing viscous friction.
##
# Factory: friction_model.

# ISA FrictionModel
from pylith.friction.FrictionModel import FrictionModel

# Import the SWIG module ViscousFriction object and rename it
# ModuleViscousFriction so that it doesn't clash with the local Python
# class of the same name.
from frictioncontrib import ViscousFriction as ModuleViscousFriction

# ViscousFriction class


class ViscousFriction(FrictionModel, ModuleViscousFriction):
    """
    Python object implementing viscous friction.

    Factory: friction_model.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="viscousfriction"):
        """
        Constructor.
        """
        FrictionModel.__init__(self, name)
        # Set the fields that are available for output. These are the
        # stored physical properties and state variables. The friction
        # model information is output with the fault information, so we
        # can also output slip, slip rate and the fault tractions.
        #
        # There are no cell fields because the fault constitutive model
        # operations on quantities evaluated at the fault vertices.
        #
        # Do not change the name of this variable. The output manager will
        # request this variable by name.
        self.availableFields = \
            {'vertex':
             {'info': ["static_coefficient",
                       "reference_slip_rate"],
              'data': ["slip_rate"]},
             'cell':
             {'info': [],
              'data': []}}
        self._loggingPrefix = "FrVisc "  # Prefix that appears in PETSc logging
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object. This
        function is called automatically by the generic Python FrictionModel
        object. It must have this name and self as the only argument.
        """
        ModuleViscousFriction.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

# This is the function that is called when you invoke
# friction = pylith.pylith.contrib.ViscousFfriction
# The name of this function MUST be 'friction_model'.
def friction_model():
    """
    Factory associated with ViscousFriction.
    """
    return ViscousFriction()  # Return our object


# End of file
