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
# Initial attempt to compute plane strain Jacobian matrices symbolically.
# PREREQUISITES:  sympy
# ----------------------------------------------------------------------
#
# import pdb
# pdb.set_trace()
import sympy
import sympy.tensor
import sympy.tensor.array
from sympy.matrices import matrix_multiply_elementwise
from itertools import product
import jacobian_common as jc
# ----------------------------------------------------------------------
print('Defining common arrays:')

# Whether to check symmetries.
checkJacobianSymmetry = True

# Constants.
zero = sympy.sympify(0)
one = sympy.sympify(1)
two = sympy.sympify(2)
three = sympy.sympify(3)
four = sympy.sympify(4)
ndim = 2

# Define basis and displacement vector.
from sympy.abc import x, y
X = [x, y]
u1 = sympy.Function('u1')(x,y)
u2 = sympy.Function('u2')(x,y)
U = [u1, u2]

# Deformation gradient, transpose, and strain tensor.
defGrad = sympy.derive_by_array(U, X)
defGradTranspose = defGrad.transpose()
strain = ((defGrad + defGradTranspose)/two).tomatrix()

# Define volumetric strain and deviatoric strain.
volStrain = sympy.tensorcontraction(strain, (0, 1))
volStrainArr = volStrain * sympy.eye(ndim)
devStrain = strain - volStrainArr/three

# ----------------------------------------------------------------------
# Power-law viscoelastic.
print('Power-law viscoelastic material:')
fileName = 'elasticity-powerlaw_iso2d.txt'
aT, n, alpha, deltaT = sympy.symbols('aT n alpha deltaT')
bulkModulus, shearModulus = sympy.symbols('bulkModulus shearModulus')
s11, s12, s22 = sympy.symbols('s11 s12 s22')
s11T, s12T, s22T = sympy.symbols('s11T s12T s22T')
s11I, s12I, s22I = sympy.symbols('s11I s12I s22I')
s11FTau, s12FTau, s22FTau = sympy.symbols('s11FTau s12FTau s22FTau')
j2FTplusDt, j2FT, j2FTau, gammaFTau = sympy.symbols('j2FTplusDt j2FT j2FTau gammaFTau')
meanStress = bulkModulus * volStrainArr
devStress = sympy.Matrix([[s11, s12], [s12, s22]])
devStressT = sympy.Matrix([[s11T, s12T], [s12T, s22T]])
devStressI = sympy.Matrix([[s11I, s12I], [s12I, s22I]])
devStressTau = alpha*devStress + (one - alpha)*devStressT
j2TplusDt = sympy.sqrt(jc.innerProd(devStress, devStress)/two)
j2T = sympy.sqrt(jc.innerProd(devStressT, devStressT)/two)
j2Tau = alpha*j2TplusDt + (one - alpha)*j2T
aE = one/(two*shearModulus)
gammaTau = aT*(j2Tau)**(n-one)
G = aE + deltaT*alpha*gammaTau
H = deltaT*(one - alpha)*gammaTau
GMat = G*sympy.ones(2,2)
HMat = H*sympy.ones(2,2)

# Derivatives.
dGDdevStress = sympy.derive_by_array(G, devStress).tomatrix()
dHDdevStress = sympy.derive_by_array(H, devStress).tomatrix()
dDevStrainDstrain = sympy.derive_by_array(devStrain, defGrad)

# Denominator.
devStressArr = sympy.Array(devStress)
devStressTArr = sympy.Array(devStressT)
denom = GMat + matrix_multiply_elementwise(devStress, dGDdevStress) + matrix_multiply_elementwise(devStressT, dHDdevStress)

denomSimp = denom.subs([(gammaTau, gammaFTau),
                        (j2Tau, j2FTau),
                        (j2TplusDt, j2FTplusDt),
                        (j2T, j2FT),
                        (devStressTau[0], s11FTau),
                        (devStressTau[1], s12FTau),
                        (devStressTau[3], s22FTau)])

jacobianDev = jc.divideTensors(dDevStrainDstrain, denomSimp)
jacobianVol = sympy.derive_by_array(meanStress, defGrad)
jacobian = jacobianDev + jacobianVol
jc.writeJacobianInfo(fileName, jacobian)
if (checkJacobianSymmetry):
    jc.checkSymmetry(jacobian, 'powerlaw2d_jacobian_symmetry.txt')
