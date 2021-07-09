#!/usr/bin/env python3
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
# Initial attempt to compute governing equations for Maxwell 3D
# material and then run assumed solution through these equations.
# PREREQUISITES:  sympy
# ----------------------------------------------------------------------
#
from sympy.abc import x, y, z, t
import sympy
import sympy.tensor
import sympy.tensor.array
# ----------------------------------------------------------------------
ndim = 3
numComps = 3
ndimRange = range(ndim)
numCompsRange = range(numComps)
outFile = 'mms_max3d.txt'

zero = sympy.sympify(0)
one = sympy.sympify(1)
two = sympy.sympify(2)
three = sympy.sympify(3)

out = open(outFile, 'w')

# ----------------------------------------------------------------------


def printTensor(tensor, tensorName):
    """Function to print components of a tensor.
    For now, assume a rank 1 or 2 tensor.
    """

    print("  Printing tensor %s" % tensorName)
    rank = tensor.rank()
    ndim = tensor.shape[0]
    simpTensor = sympy.simplify(tensor)

    for i in range(ndim):
        if (rank == 2):
            for j in range(ndim):
                line = tensorName + \
                    '_%d%d = %s\n' % (i + 1, j + 1, simpTensor[i, j])
                out.write(line)
        else:
            line = tensorName + '_%d = %s\n' % (i + 1, simpTensor[i])
            out.write(line)

    out.write('\n')

    return
# ----------------------------------------------------------------------


print("Defining basis, solution, and constants:")
# Define basis and displacement vector.
u1, u2, u3 = sympy.symbols('u1 u2 u3', type="Function")
X = sympy.tensor.array.Array([x, y, z])
dt = sympy.symbols('dt')

# Material constants.
(bulkModulus, shearModulus,
 maxwellTime) = sympy.symbols('bulkModulus shearModulus maxwellTime')

# Assumed displacements:  second-order spatial variation for now.
a, b, c, d, e, f, g = sympy.symbols('a b c d e f g')
u1 = (x * x * a + two * x * y * b + y * y * c + two * x * z * d + two * y * z * e + z * z * f) * \
    sympy.exp(-t / maxwellTime)
u2 = (x * x * c + two * x * y * b + y * y * a + two * x * z * d + two * y * z * e + z * z * f) * \
    sympy.exp(-t / maxwellTime)
u3 = (x * x * f + two * x * y * b + y * y * c + two * x * z * d + two * y * z * e + z * z * a) * \
    sympy.exp(-t / maxwellTime)
U = sympy.tensor.array.Array([u1, u2, u3])
Udot = U.diff(t)

# Perturbed solution:
u1Pert = u1 + g * x
u2Pert = u2 + g * x
u3Pert = u3 + g * x
UPert = sympy.tensor.array.Array([u1Pert, u2Pert, u3Pert])
UdotPert = UPert.diff(t)

# Deformation gradient, transpose, and strain tensor.
print("Computing strain tensor:")
defGrad = sympy.tensor.array.derive_by_array(U, X)
defGradTranspose = sympy.tensor.array.Array(defGrad.tomatrix().transpose())
strain = (defGrad + defGradTranspose) / two
strainRate = strain.diff(t)

# Deformation gradient, etc., for perturbed solution.
print("Computing strain tensor for perturbed solution:")
defGradPert = sympy.tensor.array.derive_by_array(UPert, X)
defGradTransposePert = sympy.tensor.array.Array(
    defGradPert.tomatrix().transpose())
strainPert = (defGradPert + defGradTransposePert) / two
strainRatePert = strainPert.diff(t)

# Define volumetric strain and deviatoric strain.
print("Computing deviatoric strain tensor:")
volStrain = sympy.tensor.array.tensorcontraction(strain, (0, 1))
volStrainArr = sympy.tensor.array.tensorproduct(volStrain, sympy.eye(ndim))
devStrain = strain - volStrainArr / three
devStrainRate = devStrain.diff(t)

# Define volumetric strain and deviatoric strain for perturbed solution.
print("Computing deviatoric strain tensor for perturbed solution:")
volStrainPert = sympy.tensor.array.tensorcontraction(strainPert, (0, 1))
volStrainArrPert = sympy.tensor.array.tensorproduct(
    volStrainPert, sympy.eye(ndim))
devStrainPert = strainPert - volStrainArrPert / three
devStrainRatePert = devStrainPert.diff(t)

# Assumed viscous strain.
print("Computing viscous strains:")
dq = maxwellTime * (one - sympy.exp(-t / maxwellTime)) / t
visStrain = devStrain * dq

# Assumed viscous strain for perturbed solution.
print("Computing viscous strains for perturbed solution:")
expFac = sympy.exp(-dt / maxwellTime)
dqPert = maxwellTime * (one - sympy.exp(-dt / maxwellTime)) / dt
visStrainPert = expFac * visStrain + dq * (devStrainPert - devStrain)

# Define viscous strain rate and stress function.
visStrainRate = visStrain.diff(t)
visStrainFunc = maxwellTime * (devStrainRate - visStrainRate)

# Define viscous strain rate and stress function for perturbed solution.
visStrainRatePert = visStrainPert.diff(t)
visStrainFuncPert = maxwellTime * (devStrainRatePert - visStrainRatePert)

# Define deviatoric stress and mean stress.
print("Computing stresses:")
devStress = two * shearModulus * visStrainFunc
meanStressArr = bulkModulus * volStrainArr
stress = meanStressArr + devStress

# Define deviatoric stress and mean stress for perturbed solution.
print("Computing stresses for perturbed solution:")
devStressPert = two * shearModulus * visStrainFuncPert
meanStressArrPert = bulkModulus * volStrainArrPert
stressPert = meanStressArrPert + devStressPert

# Equilibrium equation.
print("Computing equilibrium:")
equilDeriv = sympy.tensor.array.derive_by_array(stress, X)
equil = sympy.tensor.array.tensorcontraction(equilDeriv, (1, 2))

# Write results to file.
print("Writing solution variables:")
out.write('Solution variables:\n')
printTensor(U, 's')
printTensor(Udot, 's_t')
printTensor(defGrad, 's_x')

out.write('Solution variables (perturbed solution):\n')
printTensor(UPert, 'sPert')
printTensor(UdotPert, 'sPert_t')
printTensor(defGradPert, 'sPert_x')

print("Writing auxiliary variables:")
out.write('\nAuxiliary variables:\n')
printTensor(strain, 'totalStrain')
printTensor(strainRate, 'totalStrain_t')
printTensor(visStrain, 'visStrain')
printTensor(visStrainRate, 'visStrain_t')

out.write('\nAuxiliary variables (perturbed solution):\n')
printTensor(strainPert, 'totalStrainPert')
printTensor(strainRatePert, 'totalStrainPert_t')
printTensor(visStrainPert, 'visStrainPert')
printTensor(visStrainRatePert, 'visStrainPert_t')

print("Writing equilibrium equations:")
out.write('\nEquilibrium:\n')
printTensor(equil, 'equil')

print("Writing additional variables:")
out.write('\nAdditional:\n')
printTensor(devStrain, 'devStrain')
printTensor(devStrainRate, 'devStrain_t')
printTensor(stress, 'stress')
printTensor(devStress, 'devStress')
printTensor(visStrainFunc, 'visStrainFunc')

out.write('\nAdditional (perturbed solution):\n')
printTensor(devStrainPert, 'devStrainPert')
printTensor(devStrainRatePert, 'devStrainPert_t')
printTensor(stressPert, 'stressPert')
printTensor(devStressPert, 'devStressPert')
printTensor(visStrainFuncPert, 'visStrainFuncPert')

out.close()
