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
# @file pythia.pyre/meshio/DataWriter.py
#
# @brief Python abstract base class for writing finite-element data.
#
# Factory: output_data_writer

from pylith.utils.PetscComponent import PetscComponent

import os


class DataWriter(PetscComponent):
    """Python abstract base class for writing finite-element data.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="datawriter"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="datawriter")
        return

    def preinitialize(self):
        """Setup data writer.
        """
        self._createModuleObj()
        return

    @staticmethod
    def mkfilename(outputDir, simName, label, suffix):
        """Create filename from output directory, simulation name, label, and filename suffix.
        """
        filename = os.path.join(outputDir, "{}-{}.{}".format(simName, label, suffix))
        return filename

    def mkpath(self, filename):
        """Create path for output file.
        """
        self._info.log("Creating path for output file '{}'".format(filename))
        relpath = os.path.dirname(filename)

        if relpath and not os.path.exists(relpath):
            # Only create directory on proc 0
            from pylith.mpi.Communicator import mpi_comm_world
            comm = mpi_comm_world()
            if not comm.rank:
                os.makedirs(relpath)
        return

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        return

    # PRIVATE METHODS /////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Create handle to C++ object."""
        raise NotImplementedError("Implement in subclass.")


# End of file
