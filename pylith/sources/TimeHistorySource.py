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
# @file pylith/sources/TimeHistorySource.py
#
# @brief Python object for solving the timehistorysource equation.
#
# Factory: source

from .Source import Source
from .sources import TimeHistorySource as ModuleTimeHistorySource


class TimeHistorySource(Source, ModuleTimeHistorySource):
    """Python source property manager.

    FACTORY: source
    """

    import pythia.pyre.inventory

    useInitial = pythia.pyre.inventory.bool("use_initial", default=True)
    useInitial.meta['tip'] = "Use initial term in time-dependent expression."

    useRate = pythia.pyre.inventory.bool("use_rate", default=False)
    useRate.meta['tip'] = "Use rate term in time-dependent expression."

    useTimeHistory = pythia.pyre.inventory.bool("use_time_history", default=True)
    useTimeHistory.meta['tip'] = "Use time history term in time-dependent expression."

    dbTimeHistory = pythia.pyre.inventory.facility(
        "time_history", factory=NullComponent, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude as a function of time."
    
    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="timehistorysource"):
        """Constructor.
        """
        Source.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsTimeHistorySource import AuxSubfieldsTimeHistorySource
        self.auxiliarySubfields = AuxSubfieldsTimeHistorySource("auxiliary_subfields")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_is_root
        if mpi_is_root():
            self._info.log(
                "Performing minimal initialization of time-dependent Neumann boundary condition '%s'." % self.aliases[-1])
        
        
        Source.preinitialize(self, problem)
        ModuleTimeHistorySource.useTimeHistory(self, self.useTimeHistory)
        if not isinstance(self.dbTimeHistory, NullComponent):
            ModuleTimeHistorySource.setTimeHistoryDB(
                self, self.dbTimeHistory)
        return
    
    def _validate(self, context):
        if isinstance(self.inventory.dbTimeHistory, NullComponent):
            trait = self.inventory.getTrait("time_history")
            self._validationError(context, trait,
                f"Missing time history database for time history wavelet source '{self.aliases[-1]}'.")

    def _validationError(self, context, trait, msg):
        from pythia.pyre.inventory.Item import Item
        error = ValueError(msg)
        descriptor = self.getTraitDescriptor(trait.name)
        context.error(error, items=[Item(trait, descriptor)])




    def preinitialize(self, problem):
        """Setup source.
        """
        Source.preinitialize(self, problem)


        return

    def _createModuleObj(self):
        """Create handle to C++ TimeHistorySource.
        """
        ModuleTimeHistorySource.__init__(self)
        return


# Factories

def source():
    """Factory associated with TimeHistorySource.
    """
    return TimeHistorySource()


# End of file
