# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2018 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/cornercases/3d/axialdisp_rotated_soln.py
#
# @brief Analytical solution to axial extension in x-direction.
#
# Dirichlet boundary conditions in unrotated coordinates.
# boundary_xneg
#   Ux(-4000,y,z) = -b
#   Uy(-4000,y,z) = -c(y,z)
#   Uz(-4000,y,z) = -d(y,z)
# boundary_xpos:
#   Tx(+4000,0,y,z) = +5MPa

import math
import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density * p_vs**2
p_lambda = p_density * p_vp**2 - 2 * p_mu

# Rotation angles.
alpha = numpy.radians(30.0) # Rotation about x axis.
beta = numpy.radians(30.0) # Rotation about y axis.

# Rotation matrix.
ca = math.cos(alpha)
sa = math.sin(alpha)
cb = math.cos(beta)
sb = math.sin(beta)
rotMatA = numpy.array([[ ca, -sa, 0.0],
                       [ sa,  ca, 0.0],
                       [0.0, 0.0, 1.0]], dtype=numpy.float64)
rotMatB = numpy.array([[ cb, 0.0,  sb],
                       [0.0, 1.0, 0.0],
                       [-sb, 0.0,  cb]], dtype=numpy.float64)
rotMat = numpy.dot(rotMatA, rotMatB)

# Uniform stress field
sxx = +5.0e+6
syy = 0.0
szz = 0.0
sxy = 0.0
syz = 0.0
sxz = 0.0

# Uniform strain field
exx = 1.0 / (2 * p_mu) * (sxx - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
eyy = 1.0 / (2 * p_mu) * (syy - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))
ezz = 1.0 / (2 * p_mu) * (szz - p_lambda / (3 * p_lambda + 2 * p_mu) * (sxx + syy + szz))

exy = 1.0 / (2 * p_mu) * (sxy)
eyz = 1.0 / (2 * p_mu) * (syz)
exz = 1.0 / (2 * p_mu) * (sxz)

# Stress and strain fields in rotated coordinates.
stressTensor = numpy.array([[sxx, sxy, sxz],
                            [sxy, syy, syz],
                            [sxz, syz, szz]], dtype=numpy.float64)
strainTensor = numpy.array([[exx, exy, exz],
                            [exy, eyy, eyz],
                            [exz, eyz, ezz]], dtype=numpy.float64)
stressTensorRotated = numpy.dot(numpy.dot(rotMat, stressTensor), rotMat.transpose())
strainTensorRotated = numpy.dot(numpy.dot(rotMat, strainTensor), rotMat.transpose())


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to axial extension problem.
    """
    SPACE_DIM = 3
    TENSOR_SIZE = 6

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "initial_amplitude": self.displacement,
        }
        return

    def getField(self, name, pts):
        return self.fields[name](pts)

    def displacement(self, locs):
        """
        Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        locsUnrotated = numpy.dot(locs, rotMat)
        dispRotated = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp = numpy.zeros((npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[:, 0] = exx * locsUnrotated[:, 0] + exy * locsUnrotated[:, 1] + exz * locsUnrotated[:, 2]
        disp[:, 1] = exy * locsUnrotated[:, 0] + eyy * locsUnrotated[:, 1] + eyz * locsUnrotated[:, 2]
        disp[:, 2] = exz * locsUnrotated[:, 0] + eyz * locsUnrotated[:, 1] + ezz * locsUnrotated[:, 2]
        dispRotated[0,:,:] = numpy.dot(disp, rotMat.transpose())
        return dispRotated

    def density(self, locs):
        """
        Compute density field at locations.
        """
        (npts, dim) = locs.shape
        density = p_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def shear_modulus(self, locs):
        """
        Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = p_mu * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def bulk_modulus(self, locs):
        """
        Compute bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        bulk_modulus = (p_lambda + 2.0 / 3.0 * p_mu) * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return bulk_modulus

    def strain(self, locs):
        """
        Compute strain field at locations.
        """
        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[0, :, 0] = strainTensorRotated[0,0]
        strain[0, :, 1] = strainTensorRotated[1,1]
        strain[0, :, 2] = strainTensorRotated[2,2]
        strain[0, :, 3] = strainTensorRotated[0,1]
        strain[0, :, 4] = strainTensorRotated[1,2]
        strain[0, :, 5] = strainTensorRotated[0,2]
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0, :, 0] = stressTensorRotated[0,0]
        stress[0, :, 1] = stressTensorRotated[1,1]
        stress[0, :, 2] = stressTensorRotated[2,2]
        stress[0, :, 3] = stressTensorRotated[0,1]
        stress[0, :, 4] = stressTensorRotated[1,2]
        stress[0, :, 5] = stressTensorRotated[0,2]
        return stress


# End of file
