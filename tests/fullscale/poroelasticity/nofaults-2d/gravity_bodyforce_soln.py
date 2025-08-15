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
p_bulk_density = (1.0 - p_porosity) * p_solid_density + p_porosity * p_fluid_density


gacc = 9.80665  # m/s
ymax = 0.0
ymin = -8000.0  # m
fy = 5.0e3  # Pa/m


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to gravity problem."""

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
                "bc_disp_xneg": self.displacement_zero,
                "bc_disp_xpos": self.displacement_zero,
                "bc_disp_yneg": self.displacement_zero,
                "bc_press_ypos": self.pressure_zero,
            },
            "normal_dir": {
                "bc_disp_xneg": self.orientation_dir((-1, 0)),
                "bc_disp_xpos": self.orientation_dir((+1, 0)),
                "bc_disp_yneg": self.orientation_dir((0, -1)),
                "bc_press_ypos": self.orientation_dir((+1, 0)),
            },
            "tangential_dir": {
                "bc_disp_xneg": self.orientation_dir((0, -1)),
                "bc_disp_xpos": self.orientation_dir((0, +1)),
                "bc_disp_yneg": self.orientation_dir((+1, 0)),
                "bc_press_ypos": self.orientation_dir((0, +1)),
            },
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
        y = locs[:, 1]

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        uy = (
            -0.5
            / (p_lambda + 2 * p_mu)
            * ((p_bulk_density - p_alpha * p_fluid_density) * gacc - (1 - p_alpha) * fy)
            * ((ymax - ymin) ** 2 - y**2)
        )
        disp[0, :, 0] = 0.0
        disp[0, :, 1] = uy
        return disp

    def displacement_zero(self, locs):
        (npts, _) = locs.shape

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return disp

    def pressure(self, locs):
        """Compute pressure field at locations."""
        (npts, _) = locs.shape

        y = locs[:, 1]

        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        pressure[0, :, 0] = (p_fluid_density * gacc - fy) * (ymax - y)
        return pressure

    def pressure_zero(self, locs):
        (npts, _) = locs.shape

        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        return pressure

    def trace_strain(self, locs):
        """Compute trace_strain field at locations."""
        (npts, _) = locs.shape

        y = locs[:, 1]
        p_alpha = p_biot_coefficient

        ev = (
            -1.0
            / (p_lambda + 2.0 * p_mu)
            * ((p_bulk_density - p_alpha * p_fluid_density) * gacc - (1 - p_alpha) * fy)
            * (ymax - y)
        )

        trace = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        trace[0, :, 0] = ev
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

        y = locs[:, 1]
        p_alpha = p_biot_coefficient

        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        eyy = (
            -1
            / (p_lambda + 2 * p_mu)
            * ((p_bulk_density - p_alpha * p_fluid_density) * gacc - (1 - p_alpha) * fy)
            * (ymax - y)
        )
        strain[0, :, 0] = 0.0
        strain[0, :, 1] = eyy
        strain[0, :, 2] = 0.0
        strain[0, :, 3] = 0.0
        return strain

    def stress(self, locs):
        """Compute stress field at locations."""
        (npts, _) = locs.shape

        y = locs[:, 1]
        p_alpha = p_biot_coefficient
        syy = -(p_bulk_density * gacc - fy) * (ymax - y)
        sxx = szz = -p_lambda / (p_lambda + 2 * p_mu) * (p_bulk_density * gacc - fy) * (
            ymax - y
        ) - 2 * p_mu / (p_lambda + 2 * p_mu) * p_alpha * (
            p_fluid_density * gacc - fy
        ) * (
            ymax - y
        )

        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0, :, 0] = sxx
        stress[0, :, 1] = syy
        stress[0, :, 2] = szz
        stress[0, :, 3] = 0.0
        return stress

    def orientation_dir(self, vector):
        def fn_dir(locs):
            (npts, _) = locs.shape
            values = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
            for d in range(self.SPACE_DIM):
                values[:, :, d] = vector[d]
            return values

        return fn_dir


# End of file
