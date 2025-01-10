# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import os

from pylith.utils.PetscComponent import PetscComponent


class DataWriter(PetscComponent):
    """
    Abstract base class writing solution, auxiliary, and derived subfields.
    """

    def __init__(self, name="datawriter"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="datawriter")

    def preinitialize(self):
        """Setup data writer.
        """
        self._createModuleObj()

    @staticmethod
    def mkfilename(outputDir, simName, label, suffix):
        """Create filename from output directory, simulation name, label, and filename suffix.
        """
        filename = os.path.join(outputDir, "{}-{}.{}".format(simName, label, suffix))
        return filename

    def mkpath(self, filename):
        """Create path for output file.
        """
        from pylith.mpi.Communicator import mpi_is_root
        isRoot = mpi_is_root()
        if isRoot:
            self._info.log("Creating path for output file '{}'".format(filename))
        relpath = os.path.dirname(filename)

        if relpath and not os.path.exists(relpath) and isRoot:
            os.makedirs(relpath)

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """

    def _createModuleObj(self):
        """Create handle to C++ object."""
        raise NotImplementedError("Implement in subclass.")


# End of file
