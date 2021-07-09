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
# @file pylith/topology/Subfield.py
#
# @brief Python object for defining attributes of a subfield within a
# field.
#
# Factory: subfield.

from pythia.pyre.components.Component import Component


class Subfield(Component):
    """Python object for defining discretization of a subfield.

    FACTORY: subfield
    """

    import pythia.pyre.inventory

    basisOrder = pythia.pyre.inventory.int("basis_order", default=1)
    basisOrder.meta['tip'] = "Order of basis functions."

    quadOrder = pythia.pyre.inventory.int("quadrature_order", default=-1)
    quadOrder.meta['tip'] = "Order of numerical quadrature."

    dimension = pythia.pyre.inventory.int("dimension", default=-1)
    dimension.meta["tip"] = "Topological dimension associated with subfield (=-1 will use dimension of domain)."

    cellBasisStr = pythia.pyre.inventory.str("cell_basis", default="default",
                                      validator=pythia.pyre.inventory.choice(["simplex", "tensor", "default"]))
    cellBasisStr.meta['tip'] = "Type of cell basis functions (simplex, tensor, or default). Default is to use type matching cell type."

    isBasisContinuous = pythia.pyre.inventory.bool("is_basis_continous", default=True)
    isBasisContinuous.meta['tip'] = "Is basis continuous?"

    feSpaceStr = pythia.pyre.inventory.str("finite_element_space", default="polynomial",
                                    validator=pythia.pyre.inventory.choice(["polynomial", "point"]))
    feSpaceStr.meta['tip'] = "Finite-element space (polynomial or point). Point space corresponds to delta functions at quadrature points."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="subfield"):
        """Constructor.
        """
        Component.__init__(self, name, facility="subfield")

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        from .topology import FieldBase

        Component._configure(self)

        mapBasis = {
            "simplex": FieldBase.SIMPLEX_BASIS,
            "tensor": FieldBase.TENSOR_BASIS,
            "default": FieldBase.DEFAULT_BASIS,
        }
        self.cellBasis = mapBasis[self.inventory.cellBasisStr]

        mapSpace = {
            "polynomial": FieldBase.POLYNOMIAL_SPACE,
            "point": FieldBase.POINT_SPACE,
        }
        self.feSpace = mapSpace[self.inventory.feSpaceStr]
        return


# ITEM FACTORIES ///////////////////////////////////////////////////////

def subfieldFactory(name):
    """Factory for subfield items.
    """
    from pythia.pyre.inventory import facility
    return facility(name, family="subfield", factory=Subfield)


# FACTORIES ////////////////////////////////////////////////////////////

def subfield():
    """Factory associated with Subfield.
    """
    return Subfield()


# End of file
