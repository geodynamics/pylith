#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/manual/3d/powerlaw/CylinderTest.py

## @brief Python script to test power-law implementation for steady-state solution
##        of a pressurized cylinder.

import math
import numpy
import h5py
from pylith.meshio.Xdmf import Xdmf
import pdb
pdb.set_trace()

# ----------------------------------------------------------------------
# Filenames.
dispFile = 'output/cylinder_pres_powerlaw_tet-domain.h5'
stressFile = 'output/cylinder_pres_powerlaw_tet-viscomat.h5'

# Tolerances.
velAbsTol = 1.0e-5
stressAbsTol = 1.0e-5
velRelTol = 1.0e-5
stressRelTol = 1.0e-5

# Solution parameters common to both problems.
# NOTE:  A negative difference (Pb - Pa) will yield imaginary velocity values for fractional exponents.
a = 2000.0 # Inner cylinder radius.
b = 20000.0 # Outer cylinder radius.
Pa = -1.0e8 # Pressure applied at r=a.
Pb = -1.0e7 # Pressure applied at r=b.

# PyLith solution parameters.
density = 2500.0
vS = 3464.1016
vP = 6000.0
powerLawExponent = 3.5
powerLawReferenceStrainRate = 1.0e-6
powerLawReferenceStress = 1.798919e+10

# Parameters converted to FLAC/analytical form.
AT = powerLawReferenceStrainRate/(powerLawReferenceStress**powerLawExponent)
A = 2.0*AT/(math.sqrt(3.0)**(3.0 - powerLawExponent))

# ----------------------------------------------------------------------
def computeAnalyticalVel(coords):
    """
    Compute analytical velocity solution for a given set of coordinates.
    """
    r = numpy.linalg.norm(coords[:,0:2], axis=1)
    k1 = 2.0/powerLawExponent
    k3 = (3.0/4.0)**((powerLawExponent + 1.0)/2.0)
    uDot = -A*k3*((Pb - Pa)*(k1/((b/a)**k1 - 1.0)))**powerLawExponent*(b**2.0/r)

    return (uDot)


def computeAnalyticalStress(coords):
    """
    Compute analytical stress solution for a given set of coordinates.
    """
    r = numpy.linalg.norm(coords[:,0:2], axis=1)
    k1 = 2.0/powerLawExponent
    k2 = 1.0/powerLawExponent
    k3 = (3.0/4.0)**((powerLawExponent + 1.0)/2.0)
    srr = -Pb + (Pb - Pa)*((b/r)**k1 - 1.0)/((b/a)**k1 - 1.0)
    stt = -Pb - (Pb - Pa)*((k1 - 1.0)*(b/r)**k1 + 1.0)/((b/a)**k1 - 1.0)
    szz = -Pb - (Pb - Pa)*((k2 - 1.0)*(b/r)**k1 + 1.0)/((b/a)**k1 - 1.0)

    return (srr, stt, szz)


def cylToCartStress(srr, stt, coords):
    """
    Convert stresses in cylindrical coordinates to Cartesian.
    Assumption is that shear stresses are zero.
    """
    angs = numpy.arctan2(coords[:,1], coords[:,0])
    ca = numpy.cos(angs)
    sa = numpy.sin(angs)

    sxx = srr*ca*ca + stt*sa*sa
    syy = srr*sa*sa + stt*ca*ca

    return (sxx, syy)


def cylToCartDisp(ur, coords):
    """
    Convert displacements in cylindrical coordinates to Cartesian.
    Assumption is that tangential displacements are zero.
    """
    angs = numpy.arctan2(coords[:,1], coords[:,0])
    ca = numpy.cos(angs)
    sa = numpy.sin(angs)

    ux = ur*ca
    uy = ur*sa

    return (ux, uy)


# ----------------------------------------------------------------------
# Read stress info from HDF5 file.
h5Stress = h5py.File(stressFile, 'r')
coords = h5Stress['geometry/vertices'][:]
connect = numpy.array(h5Stress['topology/cells'][:], dtype=numpy.int64)
cellCoords = coords[connect, :]
cellCenters = numpy.mean(cellCoords, axis=1)
stressNum = h5Stress['cell_fields/cauchy_stress'][-1,:,:]
h5Stress.close()

# Read displacement info from HDF5 file.
h5Disp = h5py.File(dispFile, 'r')
time = h5Disp['time'][:,0,0]
dt = time[-1] - time[-2]
dispTnMinus1 = h5Disp['vertex_fields/displacement'][-2,:,:]
dispTn = h5Disp['vertex_fields/displacement'][-1,:,:]
velNum = (dispTn - dispTnMinus1)/dt
h5Disp.close()

# Compute analytical solution.
(srrAnl, sttAnl, szzAnl) = computeAnalyticalStress(cellCenters)
urDotAnl = computeAnalyticalVel(coords)
(sxxAnl, syyAnl) = cylToCartStress(srrAnl, sttAnl, cellCenters)
(uxDotAnl, uyDotAnl) = cylToCartDisp(urDotAnl, coords)

# Compute difference.
uxDiff = uxDotAnl - velNum[:,0]
uyDiff = uyDotAnl - velNum[:,1]
uzDiff = -velNum[:,2]
sxxDiff = sxxAnl - stressNum[:,0]
syyDiff = syyAnl - stressNum[:,1]
szzDiff = szzAnl - stressNum[:,2]
sxyDiff = -stressNum[:,3]
syzDiff = -stressNum[:,4]
sxzDiff = -stressNum[:,5]

# End of file 
