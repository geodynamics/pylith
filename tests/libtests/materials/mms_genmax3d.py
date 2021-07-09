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
# Initial attempt to compute governing equations for generalized Maxwell
# 3D material and then run assumed solution through these equations.
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
outFile = 'mms_genmax3d.txt'

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
(bulkModulus, shearModulus) = sympy.symbols('bulkModulus shearModulus')
(maxwellTime_1, maxwellTime_2, maxwellTime_3) = sympy.symbols(
    'maxwellTime_1 maxwellTime_2 maxwellTime_3')
(shearModulusRatio_1, shearModulusRatio_2, shearModulusRatio_3) = sympy.symbols(
    'shearModulusRatio_1 shearModulusRatio_2 shearModulusRatio_3')

# Assumed displacements:  second-order spatial variation for now.
a, b, c, d, e, f, g = sympy.symbols('a b c d e f g')
u1 = (x * x * a + two * x * y * b + y * y * c +
      two * x * z * d + two * y * z * e + z * z * f) * \
     (shearModulusRatio_1 * sympy.exp(-t / maxwellTime_1) +
      shearModulusRatio_2 * sympy.exp(-t / maxwellTime_2) +
      shearModulusRatio_3 * sympy.exp(-t / maxwellTime_3))
u2 = (x * x * c + two * x * y * b + y * y * a +
      two * x * z * d + two * y * z * e + z * z * f) * \
     (shearModulusRatio_1 * sympy.exp(-t / maxwellTime_1) +
      shearModulusRatio_2 * sympy.exp(-t / maxwellTime_2) +
      shearModulusRatio_3 * sympy.exp(-t / maxwellTime_3))
u3 = (x * x * f + two * x * y * b + y * y * c +
      two * x * z * d + two * y * z * e + z * z * a) * \
     (shearModulusRatio_1 * sympy.exp(-t / maxwellTime_1) +
      shearModulusRatio_2 * sympy.exp(-t / maxwellTime_2) +
      shearModulusRatio_3 * sympy.exp(-t / maxwellTime_3))
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
print("Computing strain tensor perturbed solution:")
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
volStrainArrPert = sympy.tensor.array.tensorproduct(volStrainPert,
                                                    sympy.eye(ndim))
devStrainPert = strainPert - volStrainArrPert / three
devStrainRatePert = devStrainPert.diff(t)

# Assumed viscous strains.
print("Computing viscous strains:")
dq_1 = maxwellTime_1 * (one - sympy.exp(-t / maxwellTime_1)) / t
dq_2 = maxwellTime_2 * (one - sympy.exp(-t / maxwellTime_2)) / t
dq_3 = maxwellTime_3 * (one - sympy.exp(-t / maxwellTime_3)) / t
visStrain_1 = devStrain * dq_1
visStrain_2 = devStrain * dq_2
visStrain_3 = devStrain * dq_3

# Assumed viscous strains for perturbed solution.
print("Computing viscous strains for perturbed solution:")
expFac_1Pert = sympy.exp(-dt / maxwellTime_1)
expFac_2Pert = sympy.exp(-dt / maxwellTime_2)
expFac_3Pert = sympy.exp(-dt / maxwellTime_3)
dq_1Pert = maxwellTime_1 * (one - sympy.exp(-dt / maxwellTime_1)) / dt
dq_2Pert = maxwellTime_2 * (one - sympy.exp(-dt / maxwellTime_2)) / dt
dq_3Pert = maxwellTime_3 * (one - sympy.exp(-dt / maxwellTime_3)) / dt
visStrain_1Pert = expFac_1Pert * visStrain_1 + \
    dq_1Pert * (devStrainPert - devStrain)
visStrain_2Pert = expFac_2Pert * visStrain_2 + \
    dq_2Pert * (devStrainPert - devStrain)
visStrain_3Pert = expFac_3Pert * visStrain_3 + \
    dq_3Pert * (devStrainPert - devStrain)

# Define viscous strain rate and stress function.
print("Computing viscous strain rates:")
visStrainRate_1 = visStrain_1.diff(t)
visStrainFunc_1 = maxwellTime_1 * (devStrainRate - visStrainRate_1)
visStrainRate_2 = visStrain_2.diff(t)
visStrainFunc_2 = maxwellTime_2 * (devStrainRate - visStrainRate_2)
visStrainRate_3 = visStrain_3.diff(t)
visStrainFunc_3 = maxwellTime_3 * (devStrainRate - visStrainRate_3)

# Define viscous strain rate and stress function for perturbed solution.
print("Computing viscous strain rates for perturbed solution:")
visStrainRate_1Pert = visStrain_1Pert.diff(t)
visStrainFunc_1Pert = maxwellTime_1 * (devStrainRatePert - visStrainRate_1Pert)
visStrainRate_2Pert = visStrain_2Pert.diff(t)
visStrainFunc_2Pert = maxwellTime_2 * (devStrainRatePert - visStrainRate_2Pert)
visStrainRate_3Pert = visStrain_3Pert.diff(t)
visStrainFunc_3Pert = maxwellTime_3 * (devStrainRatePert - visStrainRate_3Pert)

# Define deviatoric stress and mean stress.
print("Computing stresses:")
shearModulusRatio_0 = one - shearModulusRatio_1 - shearModulusRatio_2 - \
    shearModulusRatio_3

devStress = two * shearModulus * (shearModulusRatio_0 * devStrain +
                                  shearModulusRatio_1 * visStrainFunc_1 +
                                  shearModulusRatio_2 * visStrainFunc_2 +
                                  shearModulusRatio_3 * visStrainFunc_3)
meanStressArr = bulkModulus * volStrainArr
stress = meanStressArr + devStress

# Define deviatoric stress and mean stress for perturbed solution.
print("Computing stresses for perturbed solution:")
devStressPert = two * shearModulus * \
    (shearModulusRatio_0 * devStrainPert +
     shearModulusRatio_1 * visStrainFunc_1Pert +
     shearModulusRatio_2 * visStrainFunc_2Pert +
     shearModulusRatio_3 * visStrainFunc_3Pert)
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
printTensor(visStrain_1, 'visStrain_1')
printTensor(visStrainRate_1, 'visStrain_1_t')
printTensor(visStrain_2, 'visStrain_2')
printTensor(visStrainRate_2, 'visStrain_2_t')
printTensor(visStrain_3, 'visStrain_3')
printTensor(visStrainRate_3, 'visStrain_3_t')

out.write('\nAuxiliary variables (perturbed solution):\n')
printTensor(strainPert, 'totalStrainPert')
printTensor(strainRatePert, 'totalStrainPert_t')
printTensor(visStrain_1Pert, 'visStrainPert_1')
printTensor(visStrainRate_1Pert, 'visStrainPert_1_t')
printTensor(visStrain_2Pert, 'visStrainPert_2')
printTensor(visStrainRate_2Pert, 'visStrainPert_2_t')
printTensor(visStrain_3Pert, 'visStrainPert_3')
printTensor(visStrainRate_3Pert, 'visStrainPert_3_t')

print("Writing equilibrium equations:")
out.write('\nEquilibrium:\n')
printTensor(equil, 'equil')

print("Writing additional variables:")
out.write('\nAdditional:\n')
printTensor(devStrain, 'devStrain')
printTensor(devStrainRate, 'devStrain_t')
printTensor(stress, 'stress')
printTensor(devStress, 'devStress')
printTensor(visStrainFunc_1, 'visStrainFunc_1')
printTensor(visStrainFunc_2, 'visStrainFunc_2')
printTensor(visStrainFunc_3, 'visStrainFunc_3')

out.write('\nAdditional (perturbed solution):\n')
printTensor(devStrainPert, 'devStrainPert')
printTensor(devStrainRatePert, 'devStrainPert_t')
printTensor(stressPert, 'stressPert')
printTensor(devStressPert, 'devStressPert')
printTensor(visStrainFunc_1Pert, 'visStrainFuncPert_1')
printTensor(visStrainFunc_2Pert, 'visStrainFuncPert_2')
printTensor(visStrainFunc_3Pert, 'visStrainFuncPert_3')

out.close()
