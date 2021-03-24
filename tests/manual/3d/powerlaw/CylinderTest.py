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
inFile = 'output/cylinder_pres_powerlaw_hex-viscomat.h5'
outFile = 'output/cylinder_pres_powerlaw_hex-compare_analyt.h5'

# Solution parameters common to both problems.
a = 2000.0 # Inner cylinder radius.
b = 20000.0 # Outer cylinder radius.
Pa = -1.0e7 # Pressure applied at r=a.
Pb = -1.0e8 # Pressure applied at r=b.

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
def computeAnalytical(coords):
    """
    Compute analytical solution for a given set of coordinates.
    """
    r = numpy.linalg.norm(coords[:,0:2], axis=1)
    k1 = 2.0/powerLawExponent
    k2 = 1.0/powerLawExponent
    k3 = (3.0/4.0)**((powerLawExponent + 1.0)/2.0)
    srr = -Pb + (Pb - Pa)*((b/r)**k1 - 1.0)/((b/a)**k1 - 1.0)
    stt = -Pb - (Pb - Pa)*((k1 - 1.0)*(b/r)**k1 + 1.0)/((b/a)**k1 - 1.0)
    szz = -Pb - (Pb - Pa)*((k2 - 1.0)*(b/r)**k1 + 1.0)/((b/a)**k1 - 1.0)
    uDot = -A*k3*((Pb - Pa)*(k1/((b/a)**k1 - 1.0)))**powerLawExponent*(b**2.0/r)

    return (srr, stt, szz, uDot)


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
# Read coordinates and cells from HDF5 file.
h5 = h5py.File(inFile, 'r')
coords = h5['geometry/vertices'][:]
connect = numpy.array(h5['topology/cells'][:], dtype=numpy.int64)

# Read desired fields from last time step.
stressNum = h5['cell_fields/cauchy_stress'][-1,:,:]
dispNum = h5['cell_fields/displacement'][-1,:,:]

# End of file 
