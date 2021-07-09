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
# Code to compute solution and state variables for Maxwell test.
# ----------------------------------------------------------------------
#
import math
import numpy
# ----------------------------------------------------------------------
inFile = 'tri_small.vtk'
outFile = 'tri_small_vars.vtk'

# Read and write original file.
i = open(inFile, 'r')
o = open(outFile, 'w')
lines = i.readlines()
numLines = len(lines)
for lineNum in range(numLines):
    o.write(lines[lineNum])

numPoints = 12
x = numpy.zeros(numPoints, dtype=numpy.float64)
y = numpy.zeros(numPoints, dtype=numpy.float64)
z = numpy.zeros(numPoints, dtype=numpy.float64)
pointNum = 0

# Get coordinates.
for lineNum in range(5, numPoints + 5):
    line = lines[lineNum]
    lineSplit = line.split()
    x[pointNum] = float(lineSplit[0])
    y[pointNum] = float(lineSplit[1])
    z[pointNum] = float(lineSplit[2])
    pointNum += 1

# Solution constants.
density = 4000.0
vs = 5600.0
vp = math.sqrt(3.0) * vs
viscosity = 7.91700159488e+19
shearModulus = density * vs * vs
maxwellTime = viscosity / shearModulus

a = 1.0e-7
b = 2.5e-7
c = 3.0e-7
d = 9.0e-8
t = 5.0e7
dt = 5.0e7

# Solution and perturbed solution.
dispX = (a * x * x + 2.0 * b * x * y + c * y * y) * math.exp(-t / maxwellTime)
dispY = (a * y * y + 2.0 * b * x * y + c * x * x) * math.exp(-t / maxwellTime)
dispOut = numpy.column_stack((dispX, dispY, numpy.zeros_like(dispX)))
dispPertX = dispX + d * x
dispPertY = dispY + d * x
dispPertOut = numpy.column_stack(
    (dispPertX, dispPertY, numpy.zeros_like(dispPertX)))

# Total strain and total perturbed strain.
totalStrainXX = (2.0 * a * x + 2.0 * b * y) * math.exp(-t / maxwellTime)
totalStrainYY = (2.0 * b * x + 2.0 * a * y) * math.exp(-t / maxwellTime)
totalStrainXY = (b * (x + y) + c * (x + y)) * math.exp(-t / maxwellTime)
totalStrainZZ = numpy.zeros(numPoints, dtype=numpy.float64)
totalStrainPertXX = (2.0 * a * x + 2.0 * b * y) * \
    math.exp(-t / maxwellTime) + d
totalStrainPertYY = (2.0 * b * x + 2.0 * a * y) * math.exp(-t / maxwellTime)
totalStrainPertXY = (b * (x + y) + c * (x + y)) * \
    math.exp(-t / maxwellTime) + d / 2.0
totalStrainPertZZ = numpy.zeros(numPoints, dtype=numpy.float64)

# Viscous strain.
viscousStrainTXX = 2.0 * maxwellTime * (math.exp(t / maxwellTime) - 1.0) * \
    (a * (2.0 * x - y) - b * (x - 2.0 * y)) * \
    math.exp(-2.0 * t / maxwellTime) / (3.0 * t)
viscousStrainTYY = -2.0 * maxwellTime * (math.exp(t / maxwellTime) - 1.0) * \
    (a * (x - 2.0 * y) - b * (2.0 * x - y)) * \
    math.exp(-2.0 * t / maxwellTime) / (3.0 * t)
viscousStrainTZZ = -2.0 * maxwellTime * (math.exp(t / maxwellTime) - 1.0) * \
    (a * (x + y) + b * (x + y)) * math.exp(-2.0 * t / maxwellTime) / (3.0 * t)
viscousStrainTXY = maxwellTime * (math.exp(t / maxwellTime) - 1.0) * \
    (b * (x + y) + c * (x + y)) * math.exp(-2.0 * t / maxwellTime) / t

# Deviatoric strains.
meanStrainT = (totalStrainXX + totalStrainYY) / 3.0
devStrainTXX = totalStrainXX - meanStrainT
devStrainTYY = totalStrainYY - meanStrainT
devStrainTZZ = totalStrainZZ - meanStrainT
devStrainTXY = totalStrainXY
meanStrainTplusDt = (totalStrainPertXX + totalStrainPertYY) / 3.0
devStrainTplusDtXX = totalStrainPertXX - meanStrainTplusDt
devStrainTplusDtYY = totalStrainPertYY - meanStrainTplusDt
devStrainTplusDtZZ = totalStrainPertZZ - meanStrainTplusDt
devStrainTplusDtXY = totalStrainPertXY

# Time constants.
expFac = math.exp(-dt / maxwellTime)
dq = maxwellTime * (1.0 - expFac) / dt
print("expFac = %f" % expFac)
print("dq = %f" % dq)

# Perturbed viscous strains.
viscousStrainTplusDtXX = expFac * viscousStrainTXX + \
    dq * (devStrainTplusDtXX - devStrainTXX)
viscousStrainTplusDtYY = expFac * viscousStrainTYY + \
    dq * (devStrainTplusDtYY - devStrainTYY)
viscousStrainTplusDtZZ = expFac * viscousStrainTZZ + \
    dq * (devStrainTplusDtZZ - devStrainTZZ)
viscousStrainTplusDtXY = expFac * viscousStrainTXY + \
    dq * (totalStrainPertXY - totalStrainXY)

# Write results to output file.
head = 'POINT_DATA 12\n'
headscalar = 'LOOKUP_TABLE default\n'
head1 = 'VECTORS displacement double\n'
o.write(head)
o.write(head1)
numpy.savetxt(o, dispOut)
head2 = 'VECTORS displacement_pert double\n'
o.write(head2)
numpy.savetxt(o, dispPertOut)
head3 = 'SCALARS total_strain_xx double 1\n'
o.write(head3)
o.write(headscalar)
numpy.savetxt(o, totalStrainXX)
head4 = 'SCALARS total_strain_yy double 1\n'
o.write(head4)
o.write(headscalar)
numpy.savetxt(o, totalStrainYY)
head5 = 'SCALARS total_strain_xy double 1\n'
o.write(head5)
o.write(headscalar)
numpy.savetxt(o, totalStrainXY)
head6 = 'SCALARS total_strain_pert_xx double 1\n'
o.write(head6)
o.write(headscalar)
numpy.savetxt(o, totalStrainPertXX)
head7 = 'SCALARS total_strain_pert_yy double 1\n'
o.write(head7)
o.write(headscalar)
numpy.savetxt(o, totalStrainPertYY)
head8 = 'SCALARS total_strain_pert_xy double 1\n'
o.write(head8)
o.write(headscalar)
numpy.savetxt(o, totalStrainPertXY)
head9 = 'SCALARS dev_strain_t_xx double 1\n'
o.write(head9)
o.write(headscalar)
numpy.savetxt(o, devStrainTXX)
head10 = 'SCALARS dev_strain_t_yy double 1\n'
o.write(head10)
o.write(headscalar)
numpy.savetxt(o, devStrainTYY)
head10 = 'SCALARS dev_strain_t_zz double 1\n'
o.write(head10)
o.write(headscalar)
numpy.savetxt(o, devStrainTZZ)
head11 = 'SCALARS dev_strain_t_xy double 1\n'
o.write(head11)
o.write(headscalar)
numpy.savetxt(o, devStrainTXY)
head12 = 'SCALARS dev_strain_t_plus_dt_xx double 1\n'
o.write(head12)
o.write(headscalar)
numpy.savetxt(o, devStrainTplusDtXX)
head13 = 'SCALARS dev_strain_t_plus_dt_yy double 1\n'
o.write(head13)
o.write(headscalar)
numpy.savetxt(o, devStrainTplusDtYY)
head14 = 'SCALARS dev_strain_t_plus_dt_zz double 1\n'
o.write(head14)
o.write(headscalar)
numpy.savetxt(o, devStrainTplusDtZZ)
head15 = 'SCALARS dev_strain_t_plus_dt_xy double 1\n'
o.write(head15)
o.write(headscalar)
numpy.savetxt(o, devStrainTplusDtXY)
head16 = 'SCALARS viscous_strain_t_xx double 1\n'
o.write(head16)
o.write(headscalar)
numpy.savetxt(o, viscousStrainTXX)
head17 = 'SCALARS viscous_strain_t_yy double 1\n'
o.write(head17)
o.write(headscalar)
numpy.savetxt(o, viscousStrainTYY)
head17 = 'SCALARS viscous_strain_t_zz double 1\n'
o.write(head17)
o.write(headscalar)
numpy.savetxt(o, viscousStrainTZZ)
head18 = 'SCALARS viscous_strain_t_xy double 1\n'
o.write(head18)
o.write(headscalar)
numpy.savetxt(o, viscousStrainTXY)
head19 = 'SCALARS viscous_strain_t_plus_dt_xx double 1\n'
o.write(head19)
o.write(headscalar)
numpy.savetxt(o, viscousStrainTplusDtXX)
head20 = 'SCALARS viscous_strain_t_plus_dt_yy double 1\n'
o.write(head20)
o.write(headscalar)
numpy.savetxt(o, viscousStrainTplusDtYY)
head21 = 'SCALARS viscous_strain_t_plus_dt_zz double 1\n'
o.write(head21)
o.write(headscalar)
numpy.savetxt(o, viscousStrainTplusDtZZ)
head22 = 'SCALARS viscous_strain_t_plus_dt_xy double 1\n'
o.write(head22)
o.write(headscalar)
numpy.savetxt(o, viscousStrainTplusDtXY)
o.close()
