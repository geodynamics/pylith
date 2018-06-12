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

# @file pyre/meshio/OutputSoln.py
#
# @brief Python object for managing output of finite-element
# solution information.
#
# Factory: output_manager

from .OutputManager import OutputManager
from .meshio import OutputSoln as ModuleOutputSoln


class OutputSoln(OutputManager, ModuleOutputSoln):
    """
    Python object for managing output of finite-element solution
    information.

    INVENTORY

    Properties
      - None

    Facilities
      - None

    FACTORY: output_manager
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="OutputSoln"):
        """
        Constructor.
        """
        OutputManager.__init__(self, name)
        return

    def preinitialize(self):
        """
        Do
        """
        OutputManager.preinitialize(self)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        OutputManager._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to C++ object.
        """
        ModuleOutputSoln.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
    """
    Factory associated with OutputManager.
    """
    return OutputSoln()


# End of file
