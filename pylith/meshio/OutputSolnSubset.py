# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pyre/meshio/OutputSolnSubset.py
#
# @brief Python object for managing output of finite-element solution
# information over a subdomain.
#
# Factory: observer

from .OutputManager import OutputManager
from .meshio import OutputSolnSubset as ModuleOutputSolnSubset


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for group/nodeset/pset in mesh not specified.")
    return value


class OutputSolnSubset(OutputManager, ModuleOutputSolnSubset):
    """
    Python object for managing output of finite-element solution
    information over a subdomain.

    INVENTORY

    Properties
      - *label* Name identifier for subdomain.

    Facilities
      - None

    Factory: output_manager
    """

    import pyre.inventory

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for subdomain."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputsolnsubset"):
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
        ModuleOutputSolnSubset.label(self, self.label)
        return

    def verifyConfiguration(self, mesh):
        """
        Verify compatibility of configuration.
        """
        OutputManager.verifyConfiguration(self, mesh)
        ModuleOutputSolnSubset.verifyConfiguration(self, mesh)
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
        ModuleOutputSolnSubset.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """
    Factory associated with OutputManager.
    """
    return OutputSolnSubset()


# End of file
