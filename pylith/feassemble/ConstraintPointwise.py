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

# @file pylith/feassemble/ConstraintPointwise.py
##
# @brief Python abstract base class for constraints on operator
# actions with finite-elements.

from pylith.utils.PetscComponent import PetscComponent


def validateDOF(value):
    """
    Validate list of constrained degrees of freedom.
    """
    try:
        size = len(value)
        num = map(int, value)
        for v in num:
            if v < 0:
                raise ValueError
    except:
        raise ValueError, \
            "'constrained_dof' must be a zero based list of indices of degrees of " \
            "freedom at a vertex."
    return num


# ConstraintPointwise class
class ConstraintPointwise(PetscComponent):
    """
    Python abstract base class for constraints on operator
    actions with finite-elements.

    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(PetscComponent.Inventory):
        """
        Python object for managing ConstraintPointwise facilities and properties.
        """

        # @class Inventory
        # Python object for managing ConstraintPointwise facilities and properties.
        ##
        # \b Properties
        # @li \b field Field in solution to constrain.
        ##
        # \b Facilities
        # @li \b auxiliary_fields Discretization of auxiliary fields associated with material.
        # @li \b db_auxiliary_fields Database for auxiliary fields associated with material.

        import pyre.inventory

        field = pyre.inventory.str("field", default="displacement")
        field.meta['tip'] = "Field to constrain."

        constrainedDOF = pyre.inventory.list("constrained_dof", default=[], validator=validateDOF)
        constrainedDOF.meta['tip'] = "Constrained degrees of freedom (0=1st DOF, 1=2nd DOF, etc)."

        from pylith.topology.AuxSubfield import subfieldFactory
        from pylith.utils.EmptyBin import EmptyBin
        auxFields = pyre.inventory.facilityArray("auxiliary_fields", itemFactory=subfieldFactory, factory=EmptyBin)
        auxFields.meta['tip'] = "Discretization of constraint parameters."

        from spatialdata.spatialdb.SimpleDB import SimpleDB
        auxFieldsDB = pyre.inventory.facility("db_auxiliary_fields", family="spatial_database", factory=SimpleDB)
        auxFieldsDB.meta['tip'] = "Database for constraint parameters."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="constraintpointwise"):
        """
        Constructor.
        """
        return

    def preinitialize(self, mesh):
        """
        Setup constraint.
        """
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        try:
            PetscComponent._configure(self)
            self.auxFieldsDB(self.inventory.auxFieldsDB)

        except ValueError, err:
            aliases = ", ".join(self.aliases)
            raise ValueError("Error while configuring constraint (%s):\n%s" % (aliases, err.message))
        return

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError, \
            "Please implement _createModuleOb() in derived class."


# End of file
