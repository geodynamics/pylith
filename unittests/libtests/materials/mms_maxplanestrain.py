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
# Initial attempt to compute governing equations for Maxwell plane strain
# material and then run assumed solution through these equations.
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
outFile = 'mms_maxplanestrain.txt'

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
(bulkModulus, shearModulus,
 maxwellTime) = sympy.symbols('bulkModulus shearModulus maxwellTime')

# Assumed displacements:  second-order spatial variation for now.
a, b, c, d, e = sympy.symbols('a b c d e')
u1 = (x * x * a + two * x * y * b + y * y * c) * sympy.exp(-t/maxwellTime)
u2 = (x * x * c + two * y * x * b + y * y * a) * sympy.exp(-t/maxwellTime)
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

# Assumed viscous strain.
print "Computing viscous strains:"
visStrain = devStrain * (one - sympy.exp(-t/maxwellTime))
visStrain3d = devStrain3d * (one - sympy.exp(-t/maxwellTime))

# Define viscous strain rate and stress function.
visStrainRate = visStrain.diff(t)
visStrainFunc = maxwellTime * (devStrainRate - visStrainRate)
visStrainRate3d = visStrain3d.diff(t)
visStrainFunc3d = maxwellTime * (devStrainRate3d - visStrainRate3d)

# Define deviatoric stress and mean stress.
print "Computing stresses:"
devStress = two * shearModulus * visStrainFunc
meanStressArr = bulkModulus * volStrainArr
stress = meanStressArr + devStress
devStress3d = two * shearModulus * visStrainFunc3d
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
printTensor(visStrain, 'visStrain')
printTensor(visStrainRate, 'visStrain_t')

f.write('\nAuxiliary variables (3d):\n')
printTensor(strain3d, 'totalStrain3d')
printTensor(strainRate3d, 'totalStrain3d_t')
printTensor(visStrain3d, 'visStrain3d')
printTensor(visStrainRate3d, 'visStrain3d_t')

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
printTensor(visStrainFunc, 'visStrainFunc')

f.write('\nAdditional (3d):\n')
printTensor(devStrain3d, 'devStrain3d')
printTensor(devStrainRate3d, 'devStrain3d_t')
printTensor(stress3d, 'stress3d')
printTensor(devStress3d, 'devStress3d')
printTensor(visStrainFunc3d, 'visStrainFunc3d')

f.close()
