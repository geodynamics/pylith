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

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsIsotropicLinearElasticity(PetscComponent):
    """
    Auxiliary subfields associated with the isotropic linear elastic bulk rheology.

    :::{important}
    The auxiliary subfields (internal representation of material properties) do not necessarily match the values in the spatial database.
    For example, the spatial database uses density, Vp, and Vs instead of density, shear modulus, and bulk modulus because that is how they are usually characterized in seismic velocity models.
    PyLith converts the values provided by the user in a spatial database to the internal representation stored in the auxiliary field.
    :::
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_elastic.rheology.auxiliary_fields]
            shear_modulus.basis_order = 1
            bulk_modulus.basis_order = 1
            reference_stress.basis_order = 0
            reference_strain.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    shearModulus = pythia.pyre.inventory.facility("shear_modulus", family="auxiliary_subfield", factory=Subfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    bulkModulus = pythia.pyre.inventory.facility("bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    bulkModulus.meta['tip'] = "Bulk modulus subfield."

    referenceStress = pythia.pyre.inventory.facility("reference_stress", family="auxiliary_subfield", factory=Subfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pythia.pyre.inventory.facility("reference_strain", family="auxiliary_subfield", factory=Subfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."

    def __init__(self, name="auxsubfieldsisotropiclinearelasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")

    def _configure(self):
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsIsotropicLinearElasticity.
    """
    return AuxSubfieldsIsotropicLinearElasticity()


# End of file
