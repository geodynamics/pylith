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

from spatialdata.spatialdb.UniformDB import UniformDB


class ZeroDB(UniformDB):
    """
    Special case of a `UniformDB` spatial database with uniform zero initial amplitude values for degrees of freedom.

    Implements `SpatialDB`.
    """
    DOC_CONFIG = {
        "cfg": """
            # Zero displacement boundary condition constraining the y degree of freedom on the -y boundary.
            [pylithapp.problem.bc.bc_yneg]
            constrained_dof = [1]
            label = boundary_yneg
            field = displacement

            db_auxiliary_field = pylith.bc.ZeroDB
            db_auxiliary_field.description = Dirichlet displacement boundary condition on the -y boundary
            """
    }


    import pythia.pyre.inventory

    from pythia.pyre.units.length import m
    values = ["initial_amplitude", "initial_amplitude_x", "initial_amplitude_y", "initial_amplitude_z"]
    data = [0.0, 0.0, 0.0, 0.0]

    label = pythia.pyre.inventory.str("label", default="Zero initial amplitude spatial database.")
    label.meta["tip"] = "Label for ZeroDB spatial database."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="zerodb"):
        """Constructor.
        """
        UniformDB.__init__(self, name)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based on inventory.
        """
        self.inventory.values = self.values
        self.inventory.data = self.data
        UniformDB._configure(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def spatial_database():
    """Factory associated with ZeroDB.
    """
    return ZeroDB()


# End of file
