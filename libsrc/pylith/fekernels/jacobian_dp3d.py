#!/usr/bin/env python3
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
ndim = 3

# Define basis and displacement vector.
from sympy.abc import x, y, z
X = [x, y, z]
u1 = sympy.Function('u1')(x,y,z)
u2 = sympy.Function('u2')(x,y,z)
u3 = sympy.Function('u3')(x,y,z)
U = [u1, u2, u3]

# Deformation gradient, transpose, and strain tensor.
defGrad = sympy.derive_by_array(U, X)
defGradTranspose = defGrad.transpose()
strain = ((defGrad + defGradTranspose)/two).tomatrix()
strainTensor = (defGrad + defGradTranspose)/two

# Define volumetric strain and deviatoric strain.
volStrain = sympy.tensorcontraction(strain, (0, 1))
volStrainArr = volStrain * sympy.eye(ndim)
meanVolStrain = volStrain/three
meanVolStrainArr = volStrainArr/three
devStrain = strain - meanVolStrainArr

# ----------------------------------------------------------------------
# Drucker-Prager elastoplastic.
print('Drucker-Prager elastoplastic material:')
fileName = 'elasticity-drucker-prager_iso3d.txt'
fileName2 = 'elasticity-drucker-prager_iso3d_tensile.txt'
alphaFlow, alphaYield, beta, lam = sympy.symbols('alphaFlow alphaYield beta lam')
bulkModulus, shearModulus = sympy.symbols('bulkModulus shearModulus')
s11, s12, s13, s22, s23, s33 = sympy.symbols('s11 s12 s13 s22 s23 s33')
s11T, s12T, s13T, s22T, s23T, s33T = sympy.symbols('s11T s12T s13T s22T s23T s33T')
s11I, s12I, s13I, s22I, s23I, s33I = sympy.symbols('s11I s12I s13I s22I s23I s33I')
e11T, e12T, e13T, e22T, e23T, e33T = sympy.symbols('e11T e12T e13T e22T e23T e33T')
e11I, e12I, e13I, e22I, e23I, e33I = sympy.symbols('e11I e12I e13I e22I e23I e33I')
presI = sympy.symbols('presI')
volStressArrI = presI * sympy.eye(ndim)
devStress = sympy.Matrix([[s11, s12, s13],
                          [s12, s22, s23],
                          [s13, s23, s33]])
devStressT = sympy.Matrix([[s11T, s12T, s13T],
                           [s12T, s22T, s23T],
                           [s13T, s23T, s33T]])
devStressI = sympy.Matrix([[s11I, s12I, s13I],
                           [s12I, s22I, s23I],
                           [s13I, s23I, s33I]])
devStrainT = sympy.Matrix([[e11T, e12T, e13T],
                           [e12T, e22T, e23T],
                           [e13T, e23T, e33T]])
devStrainI = sympy.Matrix([[e11I, e12I, e13I],
                           [e12I, e22I, e23I],
                           [e13I, e23I, e33I]])
epp = devStrain - devStrainT - devStrainI
d = sympy.sqrt(jc.innerProd(epp, epp))
aE = one/(two*shearModulus)
aM = one/(three*bulkModulus)

# Lambda value.
print ("Defining lambda:")
lamDenom = 6*alphaYield*alphaFlow*aE + aM
lamNum = two*aE*aM*(three*alphaYield*meanVolStrain/aM + d/(sympy.sqrt(two)*aE) - beta + three*alphaYield*presI)
lam = lamNum/lamDenom
lamTensile = sympy.sqrt(two)*d

# Plastic strain increments.
print ("Defining plastic strain increments:")
deltaEpDev = lam*(epp + aE*devStressI)/(sympy.sqrt(two)*d)
deltaEpVol = lam*alphaFlow
deltaEpVolArr = deltaEpVol*sympy.eye(ndim)
deltaEpDevTensile = lamTensile*(epp + aE*devStressI)/(sympy.sqrt(two)*d)
deltaEpVolTensile = lamTensile*alphaFlow
deltaEpVolTensileArr = deltaEpVolTensile*sympy.eye(ndim)

# Stress tensor.
print ("Defining stress tensor and taking derivatives:")
stress = (epp - deltaEpDev)/aE + devStressI + (meanVolStrainArr - deltaEpVolArr)/aM + volStressArrI
jacobian1 = sympy.derive_by_array(stress, defGrad)
jacobian2 = sympy.derive_by_array(stress, defGradTranspose)
jacobian = (jacobian1 + jacobian2)/two

print ("Defining stress tensor and taking derivatives for tensile case:")
stressTensile = (epp - deltaEpDevTensile)/aE + devStressI + (meanVolStrainArr - deltaEpVolTensileArr)/aM + volStressArrI
jacobianTensile1 = sympy.derive_by_array(stressTensile, defGrad)
jacobianTensile2 = sympy.derive_by_array(stressTensile, defGradTranspose)
jacobianTensile = (jacobianTensile1 + jacobianTensile2)/two

# Substitution variables.
print ("Making initial variable substitutions:")
epp_xx, epp_yy, epp_zz, epp_xy, epp_yz, epp_xz = sympy.symbols('epp_xx epp_yy epp_zz epp_xy epp_yz epp_xz')
e_xx, e_yy, e_zz, e_xy, e_yz, e_xz = sympy.symbols('e_xx e_yy e_zz e_xy e_yz e_xz')
eps_xx, eps_yy, eps_zz, eps_xy, eps_yz, eps_xz = sympy.symbols('eps_xx eps_yy eps_zz eps_xy eps_yz eps_xz')
deltaEpDev_xx, deltaEpDev_yy, deltaEpDev_zz, deltaEpDev_xy, deltaEpDev_yz, deltaEpDev_xz = \
    sympy.symbols('deltaEpDev_xx deltaEpDev_yy deltaEpDev_zz deltaEpDev_xy deltaEpDev_yz deltaEpDev_xz')
lamS, deltaEpVolS, lamTensileS, deltaEpVolTensileS = sympy.symbols('lamS deltaEpVolS lamTensileS deltaEpVolTensileS')
L2p = jc.innerProd(epp, epp)/two
L2pS = sympy.symbols('L2pS')
eppArr = sympy.Matrix([[epp_xx, epp_xy, epp_xz],
                       [epp_xy, epp_yy, epp_yz],
                       [epp_xz, epp_yz, epp_zz]])
L2pp = jc.innerProd(eppArr, eppArr)/two
L2pp2 = jc.innerProd(eppArr, eppArr)
L2ppS, L2pp2S = sympy.symbols('L2ppS L2pp2S')

substitutions = [(lam, lamS),
                 (lamTensile, lamTensileS),
                 (deltaEpVolTensile, deltaEpVolTensileS),
                 (L2pp, L2ppS),
                 (L2pp2, L2pp2S),
                 (L2p, L2pS),
                 (deltaEpVol, deltaEpVolS),
                 (epp[0,0], epp_xx),
                 (epp[1,1], epp_yy),
                 (epp[2,2], epp_zz),
                 (epp[0,1], epp_xy),
                 (epp[0,2], epp_xz),
                 (epp[1,2], epp_yz),
                 (deltaEpDev[0,0], deltaEpDev_xx),
                 (deltaEpDev[1,1], deltaEpDev_yy),
                 (deltaEpDev[2,2], deltaEpDev_zz),
                 (deltaEpDev[0,1], deltaEpDev_xy),
                 (deltaEpDev[0,2], deltaEpDev_xz),
                 (deltaEpDev[1,2], deltaEpDev_yz),
                 (devStrain[0,0], e_xx),
                 (devStrain[1,1], e_yy),
                 (devStrain[2,2], e_zz),
                 (devStrain[0,1], e_xy),
                 (devStrain[0,2], e_xz),
                 (devStrain[1,2], e_yz),
                 (strain[0,0], eps_xx),
                 (strain[1,1], eps_yy),
                 (strain[2,2], eps_zz),
                 (strain[0,1], eps_xy),
                 (strain[0,2], eps_xz),
                 (strain[1,2], eps_yz)]
jacobianSimp = jacobian.subs(substitutions)

print ("Making initial variable substitutions for tensile case:")
jacobianTensileSimp = jacobianTensile.subs(substitutions)

print ("Making secondary variable substitutions:")
jacobianSimp2 = jacobianSimp.subs([(L2pp2, L2pp2S),
                                   (L2pp, L2ppS)])

print ("Making secondary variable substitutions for tensile case:")
jacobianTensileSimp2 = jacobianTensileSimp.subs([(L2pp2, L2pp2S),
                                                 (L2pp, L2ppS)])
                             
print ("Writing Jacobian information:")
jc.writeJacobianInfo(fileName, jacobianSimp2)

print ("Writing Jacobian information for tensile case:")
jc.writeJacobianInfo(fileName2, jacobianTensileSimp2)
if (checkJacobianSymmetry):
    print ("Writing Jacobian symmetry information:")
    jc.checkSymmetry(jacobianSimp2, 'drucker-prager_ep_jacobian_symmetry.txt')

    print ("Writing Jacobian symmetry information for tensile case:")
    jc.checkSymmetry(jacobianTensileSimp2, 'drucker-prager_ep_tensile_jacobian_symmetry.txt')
