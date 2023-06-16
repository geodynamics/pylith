# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from pylith.utils.PetscComponent import PetscComponent
from .faults import KinSrc as ModuleKinSrc


class KinSrc(PetscComponent, ModuleKinSrc):
    """
    Abstract base class for a prescribed slip source.
    """

    import pythia.pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    auxFieldDB = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=SimpleDB)
    auxFieldDB.meta['tip'] = "Database for slip time function parameters."

    from pythia.pyre.units.time import second
    originTime = pythia.pyre.inventory.dimensional("origin_time", default=0.0 * second)
    originTime.meta['tip'] = "Origin time for slip source."


    def __init__(self, name="kinsrc"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="eq_kinematic_src")
        return

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        self._createModuleObj()

        ModuleKinSrc.setIdentifier(self, self.aliases[-1])
        ModuleKinSrc.auxFieldDB(self, self.auxFieldDB)

        originTimeN = self.originTime / problem.normalizer.getTimeScale()
        ModuleKinSrc.setOriginTime(self, originTimeN)
        return

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        PetscComponent._configure(self)
        return

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Please implement _createModuleOb() in derived class.")
        return


# End of file
