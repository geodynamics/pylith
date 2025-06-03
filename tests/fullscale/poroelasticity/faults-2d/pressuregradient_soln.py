# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
# @brief Analytical solution to pressure gradient for poroelasticity with zero prescribed slip.
#

import numpy


# Physical properties
p_solid_density = 2500.0
p_fluid_density = 1000.0
p_fluid_viscosity = 1.0e-3
p_porosity = 0.02

p_shear_modulus = 3.0e10
p_drained_bulk_modulus = 8.0e10
p_fluid_bulk_modulus = 1.0e10
p_biot_coefficient = 0.2
p_isotropic_permeability = 1.0e-14


p_mu = p_shear_modulus
p_lambda = p_drained_bulk_modulus - 2.0 / 3.0 * p_shear_modulus


lx = 8000.0
p0 = 4.0e6
domain_x = 8.0e3
u0 = -213.33333333333
exx = u0 / lx
eyy = 0.0
ezz = 0.0
exy = 0.0

sxx = p_lambda * (exx + eyy + ezz) + 2 * p_mu * exx
syy = p_lambda * (exx + eyy + ezz) + 2 * p_mu * eyy
szz = p_lambda * (exx + eyy + ezz) + 2 * p_mu * ezz
sxy = 2 * p_mu * exy


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to shear problem."""

    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "trace_strain": self.trace_strain,
            "solid_density": self.solid_density,
            "fluid_density": self.fluid_density,
            "fluid_viscosity": self.fluid_viscosity,
            "shear_modulus": self.shear_modulus,
            "drained_bulk_modulus": self.drained_bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "biot_modulus": self.biot_modulus,
            "isotropic_permeability": self.isotropic_permeability,
            "porosity": self.porosity,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "initial_amplitude": {
                "bc_disp_xneg": self.displacement,
                "bc_disp_xpos": self.displacement,
                "bc_disp_yneg": self.displacement_zero,
                "bc_disp_ypos": self.displacement_zero,
                "bc_press_xneg": self.displacement,
                "bc_press_xpos": self.displacement,
            },
            "normal_dir": {
                "bc_disp_xneg": self.orientation_dir((-1, 0)),
                "bc_disp_xpos": self.orientation_dir((+1, 0)),
                "bc_disp_yneg": self.orientation_dir((0, -1)),
                "bc_disp_ypos": self.orientation_dir((0, +1)),
                "bc_press_xneg": self.orientation_dir((-1, 0)),
                "bc_press_xpos": self.orientation_dir((+1, 0)),
                "fault": self.orientation_dir((0, +1)),
            },
            "tangential_dir": {
                "bc_disp_xneg": self.orientation_dir((0, -1)),
                "bc_disp_xpos": self.orientation_dir((0, +1)),
                "bc_disp_yneg": self.orientation_dir((+1, 0)),
                "bc_disp_ypos": self.orientation_dir((-1, 0)),
                "bc_press_xneg": self.orientation_dir((0, -1)),
                "bc_press_xpos": self.orientation_dir((0, +1)),
            },
            "slip": self.slip,
            "traction_change": self.traction_change,
            "strike_dir": self.orientation_dir((-1, 0)),
        }

    def getField(self, name, mesh_entity, pts):
        if isinstance(self.fields[name], dict):
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def displacement(self, locs):
        """Compute displacement field at locations."""
        (npts, _) = locs.shape

        p_mu = p_shear_modulus
        p_lambda = p_drained_bulk_modulus - 2.0 / 3.0 * p_shear_modulus
        p_alpha = p_biot_coefficient
        x = locs[:, 0]

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        u0 = -0.5 * p_alpha * (p0 / (p_lambda + 2.0 * p_mu) * (domain_x * domain_x))
        disp[0, :, 0] = u0 * x / domain_x
        disp[0, :, 1] = 0.0
        return disp

    def displacement_zero(self, locs):
        (npts, _) = locs.shape

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return disp

    def pressure(self, locs):
        """Compute pressure field at locations."""
        (npts, _) = locs.shape

        x = locs[:, 0]

        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        pressure[0, :, 0] = p0 * (1.0 - x / domain_x)
        return pressure

    def trace_strain(self, locs):
        """Compute trace_strain field at locations."""
        (npts, _) = locs.shape

        p_alpha = p_biot_coefficient

        u0 = -0.5 * p_alpha * (p0 / (p_lambda + 2.0 * p_mu) * (domain_x * domain_x))

        trace = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        trace[0, :, 0] = u0 / domain_x
        return trace

    def solid_density(self, locs):
        """Compute solid density field at locations."""
        (npts, _) = locs.shape
        density = p_solid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def fluid_density(self, locs):
        """Compute fluid density field at locations."""
        (npts, _) = locs.shape
        density = p_fluid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def fluid_viscosity(self, locs):
        """Compute fluid viscosity field at locations."""
        (npts, _) = locs.shape
        viscosity = p_fluid_viscosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return viscosity

    def shear_modulus(self, locs):
        """Compute shear modulus field at locations."""
        (npts, _) = locs.shape
        shear_modulus = p_shear_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def drained_bulk_modulus(self, locs):
        """Compute drained bulk modulus field at locations."""
        (npts, _) = locs.shape
        bulk_modulus = p_drained_bulk_modulus * numpy.ones(
            (1, npts, 1), dtype=numpy.float64
        )
        return bulk_modulus

    def biot_coefficient(self, locs):
        """Compute Biot coefficient field at locations."""
        (npts, _) = locs.shape
        biot_coeff = p_biot_coefficient * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_coeff

    def biot_modulus(self, locs):
        """Compute Biot modulus field at locations."""
        (npts, _) = locs.shape
        p_solid_bulk_modulus = p_drained_bulk_modulus / (1.0 - p_biot_coefficient)
        p_biot_modulus = 1.0 / (
            p_porosity / p_fluid_bulk_modulus
            + (p_biot_coefficient - p_porosity) / p_solid_bulk_modulus
        )
        modulus = p_biot_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return modulus

    def isotropic_permeability(self, locs):
        """Compute permeability field at locations."""
        (npts, _) = locs.shape
        permeability = p_isotropic_permeability * numpy.ones(
            (1, npts, 1), dtype=numpy.float64
        )
        return permeability

    def porosity(self, locs):
        """Compute porosity field at locations."""
        (npts, _) = locs.shape
        value = p_porosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return value

    def strain(self, locs):
        """Compute strain field at locations."""
        (npts, _) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[0, :, 0] = exx
        strain[0, :, 1] = eyy
        strain[0, :, 2] = ezz
        strain[0, :, 3] = exy
        return strain

    def stress(self, locs):
        """Compute stress field at locations."""
        (npts, _) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0, :, 0] = sxx
        stress[0, :, 1] = syy
        stress[0, :, 2] = szz
        stress[0, :, 3] = sxy
        return stress

    def slip(self, locs):
        """Compute slip field."""
        (npts, _) = locs.shape
        slip = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return slip

    def traction_change(self, locs):
        """Compute change in traction on faults."""
        (npts, _) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[:, :, 1] = syy
        return traction

    def orientation_dir(self, vector):
        def fn_dir(locs):
            (npts, _) = locs.shape
            values = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
            for d in range(self.SPACE_DIM):
                values[:, :, d] = vector[d]
            return values

        return fn_dir


# End of file
