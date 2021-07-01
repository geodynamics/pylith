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

# @file pylith/topology/MeshImporterDist.py
##
# @brief Python implementation of importing a mesh that is already
# partitioned (distributed).
##
# Factory: mesh_generator.

from .MeshGenerator import MeshGenerator


class MeshImporterDist(MeshGenerator):
    """Python implementation of importing a mesh.

    Factory: mesh_generator.
    """

    import pythia.pyre.inventory

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    reader = pythia.pyre.inventory.facility("reader", family="mesh_io", factory=MeshIOAscii)
    reader.meta['tip'] = "Mesh reader."

    from .MeshRefiner import MeshRefiner
    refiner = pythia.pyre.inventory.facility("refiner", family="mesh_refiner", factory=MeshRefiner)
    refiner.meta['tip'] = "Mesh refiner."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="meshimporter"):
        """Constructor.
        """
        MeshGenerator.__init__(self, name)
        self._loggingPrefix = "MeIm "
        return

    def create(self, normalizer, faults=None):
        """Hook for creating mesh.
        """
        from pylith.utils.profiling import resourceUsageString

        self._setupLogging()
        logEvent = "%screate" % self._loggingPrefix
        self._eventLogger.eventBegin(logEvent)

        mesh = self.reader.read(self.debug, self.interpolate)
        if self.debug:
            mesh.view()
        self._debug.log(resourceUsageString())

        # refine mesh (if necessary)
        mesh = self.refiner.refine(mesh)

        # Nondimensionalize mesh (coordinates of vertices).
        from pylith.topology.topology import MeshOps_nondimensionalize
        MeshOps_nondimensionalize(mesh, normalizer)

        self._eventLogger.eventEnd(logEvent)
        return mesh

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_generator():
    """Factory associated with MeshImporterDist.
    """
    return MeshImporterDist()


# End of file
