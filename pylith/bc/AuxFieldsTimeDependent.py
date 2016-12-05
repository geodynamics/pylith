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

# @file pylith/materials/AuxFieldsTimeDependent.py
##
# @brief Python subfields container for isotropic, linear elasticity
# subfields.

from pylith.utils.PetscComponent import PetscComponent

# AuxFieldsTimeDependent class


class AuxFieldsTimeDependent(PetscComponent):
    """
    Python subfields container for time dependent boundary conditions.

    f(x,t) = f_0(x) + \dot{f}_1(x)(t-t_1(x)) + f_2(x)a(t-t_2(x))

    f_0(x): initial_amplitude
    \dot{f}_1(x): rate_amplitude
    t_1(x): rate_start
    f_2(x): time_history_amplitude
    t_2(x): time_history_start
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(PetscComponent.Inventory):
        """Python object for managing AuxFieldsTimeDependent
        facilities and properties.

        """

        # @class Inventory
        # Python object for managing AuxFieldsTimeDependent facilities and properties.
        ##
        # \b Properties
        # @li None
        ##
        # \b Facilities
        # @li \b initial_amplitude Initial amplitude, f_0(x), subfield.
        # @li \b rate_amplitude Rate amplitude, \dot{f}_1(x), subfield.
        # @li \b rate_start Rate start time, t_1(x), subfield.
        # @li \b time_history_amplitude Time history amplitude, f_2(x), subfield.
        # @li \b time_history_start Time history start time, t_2(s), subfield.

        import pyre.inventory

        from pylith.topology.AuxSubfield import AuxSubfield

        initialAmplitude = pyre.inventory.facility("initial_amplitude", family="auxiliary_subfield", factory=AuxSubfield)
        initialAmplitude.meta['tip'] = "Initial amplitude, f_0(x), subfield."

        rateAmplitude = pyre.inventory.facility("rate_amplitude", family="auxiliary_subfield", factory=AuxSubfield)
        rateAmplitude.meta['tip'] = "Rate amplitude, \dot{f}_1(x), subfield."

        rateStart = pyre.inventory.facility("rate_start", family="auxiliary_subfield", factory=AuxSubfield)
        rateStart.meta['tip'] = "Rate starting time, t_1(x), subfield."

        timeHistoryAmplitude = pyre.inventory.facility("time_history_amplitude", family="auxiliary_subfield", factory=AuxSubfield)
        timeHistoryAmplitude.meta['tip'] = "Time history amplitude, f_2(x). subfield"

        timeHistoryStart = pyre.inventory.facility("time_history_start", family="auxiliary_subfield", factory=AuxSubfield)
        timeHistoryStart.meta['tip'] = "Time history starting time, t_2(s), subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxfieldstimedependent"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_fields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return


# End of file
