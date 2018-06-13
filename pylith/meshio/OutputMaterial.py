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

# @file pyre/meshio/OutputMaterial.py
##
# @brief Python object for managing output of finite-element
# solution information.
##
# Factory: output_manager

from .OutputManager import OutputManager
from .meshio import OutputMaterial as ModuleOutputMaterial

# OutputMaterial class


class OutputMaterial(OutputManager, ModuleOutputMaterial):
    """
    Python object for managing output of finite-element solution
    information.

    Factory: output_manager
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputmaterial"):
        """
        Constructor.
        """
        OutputManager.__init__(self, name)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self, integrator):
        """
        Create handle to C++ object.
        """
        ModuleOutputMaterial.__init__(self, integrator)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
    """
    Factory associated with OutputManager.
    """
    return OutputMaterial()


# End of file
