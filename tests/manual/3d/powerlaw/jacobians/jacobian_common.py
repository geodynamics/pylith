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
# Common functions for computing Jacobians.
# PREREQUISITES:  sympy
# ----------------------------------------------------------------------
#
import sympy
import sympy.tensor
import sympy.tensor.array
from itertools import product
# ----------------------------------------------------------------------
def checkSymmetry(jacobian, filename):
    """
    Function to test various symmetries of Jacobian.
    """
    ndim = jacobian.shape[0]
    fmt = '%s = %s:    %s'
    fmt2 = '%s = %s:    %s\n'
    f = open(filename, 'w')
    for i in range(ndim):
        ii = i + 1
        for j in range(ndim):
            jj = j + 1
            for k in range(ndim):
                kk = k + 1
                for l in range(ndim):
                    ll = l + 1
                    # Major symmetry:
                    f1 = jacobian[k,l,i,j]
                    f1C = 'C' + str(kk) + str(ll) + str(ii) + str(jj)
                    f2 = jacobian[i,j,k,l]
                    f2C = 'C' + str(ii) + str(jj) + str(kk) + str(ll)
                    test = sympy.simplify(f1 - f2)
                    result = 'FALSE'
                    if (test == 0):
                        result = 'TRUE'
                    f.write('Major symmetry:\n')
                    f.write(fmt2 % (f1C, f2C, result))
                    print 'Major symmetry:'
                    print fmt % (f1C, f2C, result)
                    # Minor symmetry 1:
                    f1 = jacobian[j,i,k,l]
                    f1C = 'C' + str(jj) + str(ii) + str(kk) + str(ll)
                    f2 = jacobian[i,j,k,l]
                    f2C = 'C' + str(ii) + str(jj) + str(kk) + str(ll)
                    test = sympy.simplify(f1 - f2)
                    result = 'FALSE'
                    if (test == 0):
                        result = 'TRUE'
                    f.write('Minor symmetry 1:\n')
                    f.write(fmt2 % (f1C, f2C, result))
                    print 'Minor symmetry 1:'
                    print fmt % (f1C, f2C, result)
                    # Minor symmetry 2:
                    f1 = jacobian[i,j,l,k]
                    f1C = 'C' + str(ii) + str(jj) + str(ll) + str(kk)
                    f2 = jacobian[i,j,k,l]
                    f2C = 'C' + str(ii) + str(jj) + str(kk) + str(ll)
                    test = sympy.simplify(f1 - f2)
                    result = 'FALSE'
                    if (test == 0):
                        result = 'TRUE'
                    f.write('Minor symmetry 2:\n')
                    f.write(fmt2 % (f1C, f2C, result))
                    print 'Minor symmetry 2:'
                    print fmt % (f1C, f2C, result)

    f.close()
  
    return
  
def divideTensors(t1, t2):
    """
    Function to divide each element of rank 4 tensor t1 by the appropriate element
    of rank 2 tensor t2.
    """
    from sympy import MutableDenseNDimArray
    t1Rank = t1.rank()
    t1Size = t1.shape[0]
    tLen = t1Size**t1Rank
    tReturn = MutableDenseNDimArray(range(tLen), t1.shape)
    for i in range(t1Size):
        for j in range(t1Size):
            denom = t2[i,j]
            for k in range(t1Size):
                for l in range(t1Size):
                    tReturn[i,j,k,l] = t1[i,j,k,l]/denom

    return tReturn


def innerProd(mat1, mat2):
    """
    Function to compute the scalar inner product of two matrices.
    I am sure there is a much easier way to do this.
    """
    test1 = isinstance(mat1, sympy.MutableDenseMatrix)
    test2 = isinstance(mat2, sympy.MutableDenseMatrix)
    m1 = mat1.copy()
    m2 = mat2.copy()
    if (not test1):
        m1 = mat1.tomatrix()
    if (not test2):
        m2 = mat2.tomatrix()
    matMult = sympy.matrix_multiply_elementwise(m1, m2)
    rowSum = matMult * sympy.ones(matMult.shape[1], 1)
    colSum = sympy.ones(1, rowSum.shape[0]) * rowSum
    scalarProd = colSum[0]

    return scalarProd


def writeJacobianUniqueVals(f, jacobian):
    """
    Function to write unique values and assign them to variables.
    """

    # Dimension information.
    ndim = jacobian.shape[0]
    numComps = jacobian.shape[2]
    ndimRange = range(ndim)
    numCompsRange = range(numComps)

    # Unique entries in Jacobian, excluding 0.
    uniqueVals = []
    numUniqueVals = 0
    for i, j, k, l in product(numCompsRange, ndimRange, numCompsRange, ndimRange):
        testVal = jacobian[i,j,k,l]
        matchedVal = False
        if (testVal == 0):
            matchedVal = True
        else: 
            for valNum in range(numUniqueVals):
                diff = testVal - uniqueVals[valNum]
                if (diff == 0):
                    matchedVal = True
                    break
        if (matchedVal == False):
            uniqueVals.append(testVal)
            numUniqueVals += 1
                
    # uniqueVals = list(set(jacobian))
    # if (0 in uniqueVals):
        # uniqueVals.remove(0)
    # numUniqueVals = len(uniqueVals)
    uniqueValNames = numUniqueVals * [None]
    usedVals = numUniqueVals * [None]
  
    f.write("/* Unique components of Jacobian. */\n")
    outFmt = "const PylithReal %s = %s;\n"

    # Loop over Jacobian components in original PyLith order.
    ui = 0
    for i, j, k, l in product(numCompsRange, ndimRange, numCompsRange, ndimRange):
        ii = i + 1
        jj = j + 1
        kk = k + 1
        ll = l + 1
        if (jacobian[i,j,k,l] in uniqueVals and jacobian[i,j,k,l] not in usedVals):
            testInd = uniqueVals.index(jacobian[i,j,k,l])
            comp = "C" + repr(ii) + repr(jj) + repr(kk) + repr(ll)
            f.write(outFmt % (comp, jacobian[i,j,k,l]))
            uniqueValNames[testInd] = comp
            usedVals[ui] = jacobian[i,j,k,l]
            ui += 1
        if (ui == numUniqueVals):
            break

    return (uniqueVals, uniqueValNames)


def writeJacobianComments(f, jacobian):
    """
    Function to write correspondence between PETSc and PyLith Jacobian values.
    """

    f.write("/* j(f,g,df,dg) = C(f,df,g,dg)\n\n")
    outFmt = "%d:  %s = %s = %s\n"

    # Dimension information.
    ndim = jacobian.shape[0]
    numComps = jacobian.shape[2]
    ndimRange = range(ndim)
    numCompsRange = range(numComps)
  
    # Loop over Jacobian components in new order.
    ui = 0
    for i, k, j, l in product(numCompsRange, numCompsRange, ndimRange, ndimRange):
        ii = i + 1
        jj = j + 1
        kk = k + 1
        ll = l + 1
        pyComp = "C" + repr(ii) + repr(jj) + repr(kk) + repr(ll)
        peComp = "j" + repr(i) + repr(k) + repr(j) + repr(l)
        f.write(outFmt % (ui, peComp, pyComp, jacobian[i,j,k,l]))
        ui += 1
            
    f.write("*/\n\n")

    return


def writeJacobianNonzero(f, jacobian, uniqueVals, uniqueValNames):
    """
    Function to write nonzero Jacobian entries using predefined value names.
    """

    f.write("/* Nonzero Jacobian entries. */\n")
  
    outFmt = "Jg3[%d] -=  %s; /* %s */\n"

    # Dimension information.
    ndim = jacobian.shape[0]
    numComps = jacobian.shape[2]
    ndimRange = range(ndim)
    numCompsRange = range(numComps)
  
    # Loop over Jacobian components in new order.
    ui = 0
    for i, k, j, l in product(numCompsRange, numCompsRange, ndimRange, ndimRange):
        ii = i + 1
        jj = j + 1
        kk = k + 1
        ll = l + 1
        peComp = "j" + repr(i) + repr(k) + repr(j) + repr(l)
        if (jacobian[i,j,k,l] != 0):
            ind = uniqueVals.index(jacobian[i,j,k,l])
            f.write(outFmt % (ui, uniqueValNames[ind], peComp))

        ui += 1

    return


def writeJacobianInfo(fileName, jacobian):
    """
    Function to write info about Jacobian.
    """
    f = open(fileName, 'w')

    (uniqueVals, uniqueValNames) = writeJacobianUniqueVals(f, jacobian)
    writeJacobianComments(f, jacobian)
    writeJacobianNonzero(f, jacobian, uniqueVals, uniqueValNames)
    f.close()

    return
# ----------------------------------------------------------------------
                  
