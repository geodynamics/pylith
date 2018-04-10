#!/usr/bin/env python
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# Initial attempt to compute governing equations for generalized Maxwell
# plane strain material and then run assumed solution through these equations.
# PREREQUISITES:  sympy
# ----------------------------------------------------------------------
#
# import pdb
# pdb.set_trace()
import sympy
import sympy.tensor
import sympy.tensor.array
# ----------------------------------------------------------------------
ndim = 2
numComps = 2
ndimRange = range(ndim)
numCompsRange = range(numComps)
outFile = 'mms_genmaxplanestrain.txt'

zero = sympy.sympify(0)
one = sympy.sympify(1)
two = sympy.sympify(2)
three = sympy.sympify(3)

f = open(outFile, 'w')

# ----------------------------------------------------------------------
def printTensor(tensor, tensorName):
  """
  Function to print components of a tensor.
  For now, assume a rank 1 or 2 tensor.
  """
  
  print "  Printing tensor %s" % tensorName
  rank = tensor.rank()
  ndim = tensor.shape[0]
  simpTensor = sympy.simplify(tensor)
  
  for i in range(ndim):
    if (rank == 2):
      for j in range(ndim):
        line = tensorName + '_%d%d = %s\n' % (i+1, j+1, simpTensor[i,j])
        f.write(line)
    else:
      line = tensorName + '_%d = %s\n' % (i+1, simpTensor[i])
      f.write(line)

  f.write('\n')

  return
# ----------------------------------------------------------------------
  
print "Defining basis, solution, and constants:"
# Define basis and displacement vector.
from sympy.abc import x, y, z, t
u1, u2, u3 = sympy.symbols('u1 u2 u3', type="Function")
X = sympy.tensor.array.Array([x, y])
X3d = sympy.tensor.array.Array([x, y, z])

# Material constants.
(bulkModulus, shearModulus) = sympy.symbols('bulkModulus shearModulus')
(maxwellTime_1, maxwellTime_2, maxwellTime_3) = sympy.symbols(
  'maxwellTime_1 maxwellTime_2 maxwellTime_3')
(shearModulusRatio_1, shearModulusRatio_2, shearModulusRatio_3) = sympy.symbols(
  'shearModulusRatio_1 shearModulusRatio_2 shearModulusRatio_3')

# Assumed displacements:  second-order spatial variation for now.
a, b, c, d, e = sympy.symbols('a b c d e')
u1 = (x * x * a + two * x * y * b + y * y * c) * \
     (shearModulusRatio_1 * sympy.exp(-t/maxwellTime_1) + \
      shearModulusRatio_2 * sympy.exp(-t/maxwellTime_2) + \
      shearModulusRatio_3 * sympy.exp(-t/maxwellTime_3))
u2 = (x * x * c + two * y * x * b + y * y * a) * \
     (shearModulusRatio_1 * sympy.exp(-t/maxwellTime_1) + \
      shearModulusRatio_2 * sympy.exp(-t/maxwellTime_2) + \
      shearModulusRatio_3 * sympy.exp(-t/maxwellTime_3))
u3 = zero
U = sympy.tensor.array.Array([u1, u2])
Udot = U.diff(t)

U3d = sympy.tensor.array.Array([u1, u2, u3])
U3ddot = U3d.diff(t)

# Deformation gradient, transpose, and strain tensor.
print "Computing strain tensor:"
defGrad = sympy.tensor.array.derive_by_array(U, X)
defGradTranspose = sympy.tensor.array.Array(defGrad.tomatrix().transpose())
strain = (defGrad + defGradTranspose)/two
strainRate = strain.diff(t)

defGrad3d = sympy.tensor.array.derive_by_array(U3d, X3d)
defGradTranspose3d = sympy.tensor.array.Array(defGrad3d.tomatrix().transpose())
strain3d = (defGrad3d + defGradTranspose3d)/two
strainRate3d = strain3d.diff(t)

# Define volumetric strain and deviatoric strain.
print "Computing deviatoric strain tensor:"
volStrain = sympy.tensor.array.tensorcontraction(strain, (0, 1))
volStrainArr = sympy.tensor.array.tensorproduct(volStrain, sympy.eye(ndim))
devStrain = strain - volStrainArr/three
devStrainRate = devStrain.diff(t)

volStrain3d = sympy.tensor.array.tensorcontraction(strain3d, (0, 1))
volStrainArr3d = sympy.tensor.array.tensorproduct(volStrain3d, sympy.eye(3))
devStrain3d = strain3d - volStrainArr3d/three
devStrainRate3d = devStrain3d.diff(t)

# Assumed viscous strains.
print "Computing viscous strains:"
visStrain_1 = devStrain * (one - sympy.exp(-t/maxwellTime_1))
visStrain_2 = devStrain * (one - sympy.exp(-t/maxwellTime_2))
visStrain_3 = devStrain * (one - sympy.exp(-t/maxwellTime_3))

visStrain3d_1 = devStrain3d * (one - sympy.exp(-t/maxwellTime_1))
visStrain3d_2 = devStrain3d * (one - sympy.exp(-t/maxwellTime_2))
visStrain3d_3 = devStrain3d * (one - sympy.exp(-t/maxwellTime_3))

# Define viscous strain rate and stress function.
visStrainRate_1 = visStrain_1.diff(t)
visStrainFunc_1 = maxwellTime_1 * (devStrainRate - visStrainRate_1)
visStrainRate_2 = visStrain_2.diff(t)
visStrainFunc_2 = maxwellTime_2 * (devStrainRate - visStrainRate_2)
visStrainRate_3 = visStrain_3.diff(t)
visStrainFunc_3 = maxwellTime_3 * (devStrainRate - visStrainRate_3)

visStrainRate3d_1 = visStrain3d_1.diff(t)
visStrainFunc3d_1 = maxwellTime_1 * (devStrainRate3d - visStrainRate3d_1)
visStrainRate3d_2 = visStrain3d_2.diff(t)
visStrainFunc3d_2 = maxwellTime_2 * (devStrainRate3d - visStrainRate3d_2)
visStrainRate3d_3 = visStrain3d_3.diff(t)
visStrainFunc3d_3 = maxwellTime_3 * (devStrainRate3d - visStrainRate3d_3)

# Define deviatoric stress and mean stress.
print "Computing stresses:"
shearModulusRatio_0 = one - shearModulusRatio_1 - shearModulusRatio_2 - \
                      shearModulusRatio_3

devStress = two * shearModulus * (shearModulusRatio_0 * devStrain + \
                                  shearModulusRatio_1 * visStrainFunc_1 + \
                                  shearModulusRatio_2 * visStrainFunc_2 + \
                                  shearModulusRatio_3 * visStrainFunc_3)
meanStressArr = bulkModulus * volStrainArr
stress = meanStressArr + devStress

devStress3d = two * shearModulus * (shearModulusRatio_0 * devStrain3d + \
                                    shearModulusRatio_1 * visStrainFunc3d_1 + \
                                    shearModulusRatio_2 * visStrainFunc3d_2 + \
                                    shearModulusRatio_3 * visStrainFunc3d_3)
meanStressArr3d = bulkModulus * volStrainArr3d
stress3d = meanStressArr3d + devStress3d

# Equilibrium equation.
print "Computing equilibrium:"
equilDeriv = sympy.tensor.array.derive_by_array(stress, X)
equil = sympy.tensor.array.tensorcontraction(equilDeriv, (1,2))

equilDeriv3d = sympy.tensor.array.derive_by_array(stress3d, X3d)
equil3d = sympy.tensor.array.tensorcontraction(equilDeriv3d, (1,2))

# Write results to file.
print "Writing solution variables:"
f.write('Solution variables:\n')
printTensor(U, 's')
printTensor(Udot, 's_t')
printTensor(defGrad, 's_x')

f.write('Solution variables (3d):\n')
printTensor(U3d, 's3d')
printTensor(U3ddot, 's3d_t')
printTensor(defGrad3d, 's3d_x')

print "Writing auxiliary variables:"
f.write('\nAuxiliary variables:\n')
printTensor(strain, 'totalStrain')
printTensor(strainRate, 'totalStrain_t')
printTensor(visStrain_1, 'visStrain_1')
printTensor(visStrainRate_1, 'visStrain_1_t')
printTensor(visStrain_2, 'visStrain_2')
printTensor(visStrainRate_2, 'visStrain_2_t')
printTensor(visStrain_3, 'visStrain_3')
printTensor(visStrainRate_3, 'visStrain_3_t')

f.write('\nAuxiliary variables (3d):\n')
printTensor(strain3d, 'totalStrain3d')
printTensor(strainRate3d, 'totalStrain3d_t')
printTensor(visStrain3d_1, 'visStrain3d_1')
printTensor(visStrainRate3d_1, 'visStrain3d_1_t')
printTensor(visStrain3d_2, 'visStrain3d_2')
printTensor(visStrainRate3d_2, 'visStrain3d_2_t')
printTensor(visStrain3d_3, 'visStrain3d_3')
printTensor(visStrainRate3d_3, 'visStrain3d_3_t')

print "Writing equilibrium equations:"
f.write('\nEquilibrium:\n')
printTensor(equil, 'equil')

f.write('\nEquilibrium (3d):\n')
printTensor(equil3d, 'equil3d')

print "Writing additional variables:"
f.write('\nAdditional:\n')
printTensor(devStrain, 'devStrain')
printTensor(devStrainRate, 'devStrain_t')
printTensor(stress, 'stress')
printTensor(devStress, 'devStress')
printTensor(visStrainFunc_1, 'visStrainFunc_1')
printTensor(visStrainFunc_2, 'visStrainFunc_2')
printTensor(visStrainFunc_3, 'visStrainFunc_3')

f.write('\nAdditional (3d):\n')
printTensor(devStrain3d, 'devStrain3d')
printTensor(devStrainRate3d, 'devStrain3d_t')
printTensor(stress3d, 'stress3d')
printTensor(devStress3d, 'devStress3d')
printTensor(visStrainFunc3d_1, 'visStrainFunc3d_1')
printTensor(visStrainFunc3d_2, 'visStrainFunc3d_2')
printTensor(visStrainFunc3d_3, 'visStrainFunc3d_3')

f.close()
