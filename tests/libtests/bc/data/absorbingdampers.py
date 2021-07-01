#!/usr/bin/env nemesis
#
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

# @file tests/libtests/bc/data/absorbingdampers.py

# @brief Python routines for calculating values for test of absorbing
# dampers.

import math

# ----------------------------------------------------------------------


def calcTri3():
    """Calculate damping constants, residual, and Jacobian values for
    absorbing dampers on mesh with tri3 cells.
    """

    dt = 0.25
    density = 2500.0
    vp = 5000.0
    vs = 3000.0
    edgeLen = 2.0**0.5
    N0 = 0.5
    N1 = 0.5
    velX = 0.5*(N0+N1)*(N0*(1.3+1.5-1.1)+N1*(1.7+2.1-1.3)) / (2.0*dt)
    velY = 0.5*(N0+N1)*(N0*(2.1+2.4-1.8)+N1*(2.3+2.4-2.2)) / (2.0*dt)
    normal = [0.5**0.5, -0.5**0.5]

    constNormal = density*vp
    constTangential = density*vs
    dampingConsts = [abs(constNormal*normal[0] - constTangential*normal[1]),
                     abs(constNormal*normal[1] + constTangential*normal[0])]
    residualX = -dampingConsts[0]*velX*edgeLen
    residualY = -dampingConsts[1]*velY*edgeLen
    residual = [residualX, residualY,
                residualX, residualY]

    j00 = edgeLen*N0**2 / (2.0*dt)
    j01 = edgeLen*N0*N1 / (2.0*dt)
    j10 = j01
    j11 = edgeLen*N1**2 / (2.0*dt)
    jacobian = [dampingConsts[0]*j00, dampingConsts[1]*j00,
                dampingConsts[0]*j01, dampingConsts[1]*j01,
                dampingConsts[0]*j10, dampingConsts[1]*j10,
                dampingConsts[0]*j11, dampingConsts[1]*j11]

    print "Absorbing boundary for tri3 mesh"
    print "damping constants:"
    for v in dampingConsts:
        print "  %16.8e" % v
    print "values for residual:"
    for v in residual:
        print "  %16.8e" % v
    print "values for jacobian:"
    for j in jacobian:
        print "  %16.8e" % j


# ----------------------------------------------------------------------
def calcQuad4():
    """Calculate damping constants, residual, and Jacobian values for
    absorbing dampers on mesh with quad4 cells.
    """

    density = 2500.0
    vp = 5000.0
    vs = 3000.0
    dt = 0.25
    N0 = 0.5
    N1 = 0.5
    vel1X = 0.5*(N0+N1)*(N0*(1.3+1.5-1.1) + N1*(1.1+1.2-1.0)) / (2.0*dt)
    vel1Y = 0.5*(N0+N1)*(N0*(2.1+2.4-1.8) + N1*(2.0+1.6-2.4)) / (2.0*dt)
    normal1 = [-1.0, 0.0]

    vel2X = 0.5*(N0+N1)*(N0*(1.9+2.4-1.4) + N1*(2.1+2.7-1.5)) / (2.0*dt)
    vel2Y = 0.5*(N0+N1)*(N0*(2.4+2.4-2.4) + N1*(2.5+3.4-1.6)) / (2.0*dt)
    normal2 = [1.0, 0.0]

    edgeLen = 2.0

    constNormal = density*vp
    constTangential = density*vs
    dampingConsts = [abs(constNormal*normal1[0] - constTangential*normal1[1]),
                     abs(constNormal*normal1[1] + constTangential*normal1[0]),
                     abs(constNormal*normal2[0] - constTangential*normal2[1]),
                     abs(constNormal*normal2[1] + constTangential*normal2[0])]
    residual1X = -dampingConsts[0]*vel1X*edgeLen
    residual1Y = -dampingConsts[1]*vel1Y*edgeLen
    residual2X = -dampingConsts[2]*vel2X*edgeLen
    residual2Y = -dampingConsts[3]*vel2Y*edgeLen
    residual = [residual1X, residual1Y,
                residual1X, residual1Y,
                residual2X, residual2Y,
                residual2X, residual2Y]

    j00 = edgeLen*N0**2 / (2.0*dt)
    j01 = edgeLen*N0*N1 / (2.0*dt)
    j10 = j01
    j11 = edgeLen*N1**2 / (2.0*dt)
    jacobian = [dampingConsts[0]*j00, dampingConsts[1]*j00,
                dampingConsts[0]*j01, dampingConsts[1]*j01,
                dampingConsts[0]*j10, dampingConsts[1]*j10,
                dampingConsts[0]*j11, dampingConsts[1]*j11,
                dampingConsts[2]*j00, dampingConsts[3]*j00,
                dampingConsts[2]*j01, dampingConsts[3]*j01,
                dampingConsts[2]*j10, dampingConsts[3]*j10,
                dampingConsts[2]*j11, dampingConsts[3]*j11]
    print "Absorbing boundary for quad4mesh"
    print "damping constants:"
    for v in dampingConsts:
        print "  %16.8e" % v
    print "vel:"
    print "  vel1: ", vel1X, "  ", vel1Y
    print "  vel2: ", vel2X, "  ", vel2Y
    print "values for residual:"
    for v in residual:
        print "  %16.8e" % v
    print "values for jacobian:"
    for j in jacobian:
        print "  %16.8e" % j


# ----------------------------------------------------------------------
def calcTet4():
    """Calculate damping constants, residual, and Jacobian values for
    absorbing dampers on mesh with tet4 cells.
    """

    dt = 0.25
    density = 2500.0
    vp = 5000.0
    vs = 3000.0
    area = 0.5
    N0 = 1.0/3.0
    N1 = 1.0/3.0
    N2 = 1.0/3.0
    velX = (N0+N1+N2)/3.0 * (N0*(1.7+2.1-1.3)+N1 *
                             (1.5+1.8-1.2)+N2*(1.1+1.2-1.0)) / (2.0*dt)
    velY = (N0+N1+N2)/3.0 * (N0*(2.3+2.4-2.2)+N1 *
                             (2.2+2.0-2.4)+N2*(2.0+1.6-2.4)) / (2.0*dt)
    velZ = (N0+N1+N2)/3.0 * (N0*(4.4+5.2-3.6)+N1 *
                             (4.0+4.6-3.4)+N2*(3.2+3.4-3.0)) / (2.0*dt)
    normal = [-1.0, 0.0, 0.0]
    tangent1 = [0.0, -1.0, 0.0]
    tangent2 = [0.0, 0.0, 1.0]

    constNormal = density*vp
    constTangential = density*vs
    dampingConsts = [abs(constNormal*normal[0] +
                         constTangential*tangent1[0] +
                         constTangential*tangent2[0]),
                     abs(constNormal*normal[1] +
                         constTangential*tangent1[1] +
                         constTangential*tangent2[1]),
                     abs(constNormal*normal[2] +
                         constTangential*tangent1[2] +
                         constTangential*tangent2[2])]
    residualX = -dampingConsts[0]*velX*area
    residualY = -dampingConsts[1]*velY*area
    residualZ = -dampingConsts[2]*velZ*area
    residual = [residualX, residualY, residualZ,
                residualX, residualY, residualZ,
                residualX, residualY, residualZ]

    j00 = area*N0**2 / (2.0*dt)
    j01 = area*N0*N1 / (2.0*dt)
    j02 = area*N0*N2 / (2.0*dt)
    j10 = j01
    j11 = area*N1**2 / (2.0*dt)
    j12 = area*N1*N2 / (2.0*dt)
    j20 = j02
    j21 = j12
    j22 = area*N2**2 / (2.0*dt)
    jacobian = [dampingConsts[0]*j00, dampingConsts[1]*j00, dampingConsts[2]*j00,
                dampingConsts[0]*j01, dampingConsts[1] *
                j01, dampingConsts[2]*j01,
                dampingConsts[0]*j10, dampingConsts[1] *
                j10, dampingConsts[2]*j10,
                dampingConsts[0]*j11, dampingConsts[1] *
                j11, dampingConsts[2]*j11,
                dampingConsts[0]*j12, dampingConsts[1] *
                j12, dampingConsts[2]*j12,
                dampingConsts[0]*j20, dampingConsts[1] *
                j20, dampingConsts[2]*j20,
                dampingConsts[0]*j21, dampingConsts[1] *
                j21, dampingConsts[2]*j21,
                dampingConsts[0]*j22, dampingConsts[1]*j22, dampingConsts[2]*j22]

    print "Absorbing boundary for hex8 mesh"
    print "damping constants:"
    for v in dampingConsts:
        print "  %16.8e" % v
    print "values for residual:"
    for v in residual:
        print "  %16.8e" % v
    print "values for jacobian:"
    for j in jacobian:
        print "  %16.8e" % j


# ----------------------------------------------------------------------
def calcHex8():
    """Calculate damping constants, residual, and Jacobian values for
    absorbing dampers on mesh with hex8 cells.
    """

    import numpy

    dt = 0.25
    density = 2500.0
    vp = 5000.0
    vs = 3000.0
    area = 1.0
    jacobianDet = 0.5
    basis = [[0.62200847,  0.16666667,  0.0446582,   0.16666667],
             [0.16666667,  0.62200847,  0.16666667,  0.0446582],
             [0.0446582,   0.16666667,  0.62200847,  0.16666667],
             [0.16666667,  0.0446582,   0.16666667,  0.62200847]]
    cells = [[2, 8, 6, 0],
             [4, 10, 8, 2]]
    dispTmdt = [[1.0,  2.4,  3.0],
                [1.1,  2.2,  3.2],
                [1.2,  2.0,  3.4],
                [1.3,  1.8,  3.6],
                [1.4,  1.6,  3.8],
                [1.5,  1.4,  4.0],
                [1.6,  1.2,  4.2],
                [1.7,  1.0,  4.4],
                [1.8,  0.8,  4.6],
                [1.9,  0.6,  4.8],
                [2.0,  0.4,  5.0],
                [2.1,  0.2,  5.2]]
    dispT = [[1.1,  2.3,  3.2],
             [1.3,  2.1,  3.6],
             [1.5,  1.9,  4.0],
             [1.7,  1.7,  4.4],
             [1.9,  1.5,  4.8],
             [2.1,  1.3,  5.2],
             [2.3,  1.1,  5.6],
             [2.5,  0.9,  6.0],
             [2.7,  0.7,  6.4],
             [2.9,  0.5,  6.8],
             [3.1,  0.3,  7.2],
             [3.3,  0.1,  7.6]]
    dispIncr = [[1.2,  1.1,  3.4],
                [1.5,  1.0,  4.0],
                [1.8,  0.9,  4.6],
                [2.1,  0.8,  5.2],
                [2.4,  0.7,  5.8],
                [2.7,  0.6,  6.4],
                [3.0,  0.5,  7.0],
                [3.3,  0.4,  7.6],
                [3.6,  0.3,  8.2],
                [3.9,  0.2,  8.8],
                [4.2,  0.1,  9.4],
                [4.5,  0.0, 10.0]]
    normal = [0.0, 1.0, 0.0]
    tangent1 = [-1.0, 0.0, 0.0]
    tangent2 = [0.0, 0.0, 1.0]

    constNormal = density*vp
    constTangential = density*vs
    dampingConsts = [abs(constNormal*normal[0] +
                         constTangential*tangent1[0] +
                         constTangential*tangent2[0]),
                     abs(constNormal*normal[1] +
                         constTangential*tangent1[1] +
                         constTangential*tangent2[1]),
                     abs(constNormal*normal[2] +
                         constTangential*tangent1[2] +
                         constTangential*tangent2[2])]

    residual = numpy.zeros((12, 3), dtype=numpy.float64)
    jacobian = numpy.zeros((12, 3, 12, 3), dtype=numpy.float64)
    for cell in cells:
        for b in basis:
            N0 = b[0]
            N1 = b[1]
            N2 = b[2]
            N3 = b[3]
            d0x = (dispT[cell[0]][0] + dispIncr[cell[0]]
                   [0] - dispTmdt[cell[0]][0])/(2.0*dt)
            d0y = (dispT[cell[0]][1] + dispIncr[cell[0]]
                   [1] - dispTmdt[cell[0]][1])/(2.0*dt)
            d0z = (dispT[cell[0]][2] + dispIncr[cell[0]]
                   [2] - dispTmdt[cell[0]][2])/(2.0*dt)
            d1x = (dispT[cell[1]][0] + dispIncr[cell[1]]
                   [0] - dispTmdt[cell[1]][0])/(2.0*dt)
            d1y = (dispT[cell[1]][1] + dispIncr[cell[1]]
                   [1] - dispTmdt[cell[1]][1])/(2.0*dt)
            d1z = (dispT[cell[1]][2] + dispIncr[cell[1]]
                   [2] - dispTmdt[cell[1]][2])/(2.0*dt)
            d2x = (dispT[cell[2]][0] + dispIncr[cell[2]]
                   [0] - dispTmdt[cell[2]][0])/(2.0*dt)
            d2y = (dispT[cell[2]][1] + dispIncr[cell[2]]
                   [1] - dispTmdt[cell[2]][1])/(2.0*dt)
            d2z = (dispT[cell[2]][2] + dispIncr[cell[2]]
                   [2] - dispTmdt[cell[2]][2])/(2.0*dt)
            d3x = (dispT[cell[3]][0] + dispIncr[cell[3]]
                   [0] - dispTmdt[cell[3]][0])/(2.0*dt)
            d3y = (dispT[cell[3]][1] + dispIncr[cell[3]]
                   [1] - dispTmdt[cell[3]][1])/(2.0*dt)
            d3z = (dispT[cell[3]][2] + dispIncr[cell[3]]
                   [2] - dispTmdt[cell[3]][2])/(2.0*dt)
            velX = N0*d0x + N1*d1x + N2*d2x + N3*d3x
            velY = N0*d0y + N1*d1y + N2*d2y + N3*d3y
            velZ = N0*d0z + N1*d1z + N2*d2z + N3*d3z

            residualX = -dampingConsts[0] * velX * area * jacobianDet
            residualY = -dampingConsts[1] * velY * area * jacobianDet
            residualZ = -dampingConsts[2] * velZ * area * jacobianDet
            residual[cell[0], :] += N0 * \
                numpy.array([residualX, residualY, residualZ])
            residual[cell[1], :] += N1 * \
                numpy.array([residualX, residualY, residualZ])
            residual[cell[2], :] += N2 * \
                numpy.array([residualX, residualY, residualZ])
            residual[cell[3], :] += N3 * \
                numpy.array([residualX, residualY, residualZ])

        for b in basis:
            j00 = jacobianDet*area*b[0]*b[0] / (2.0*dt)
            j01 = jacobianDet*area*b[0]*b[1] / (2.0*dt)
            j02 = jacobianDet*area*b[0]*b[2] / (2.0*dt)
            j03 = jacobianDet*area*b[0]*b[3] / (2.0*dt)

            j10 = jacobianDet*area*b[1]*b[0] / (2.0*dt)
            j11 = jacobianDet*area*b[1]*b[1] / (2.0*dt)
            j12 = jacobianDet*area*b[1]*b[2] / (2.0*dt)
            j13 = jacobianDet*area*b[1]*b[3] / (2.0*dt)

            j20 = jacobianDet*area*b[2]*b[0] / (2.0*dt)
            j21 = jacobianDet*area*b[2]*b[1] / (2.0*dt)
            j22 = jacobianDet*area*b[2]*b[2] / (2.0*dt)
            j23 = jacobianDet*area*b[2]*b[3] / (2.0*dt)

            j30 = jacobianDet*area*b[3]*b[0] / (2.0*dt)
            j31 = jacobianDet*area*b[3]*b[1] / (2.0*dt)
            j32 = jacobianDet*area*b[3]*b[2] / (2.0*dt)
            j33 = jacobianDet*area*b[3]*b[3] / (2.0*dt)

            jj = [dampingConsts[0]*j00, dampingConsts[1]*j00, dampingConsts[2]*j00,
                  dampingConsts[0]*j01, dampingConsts[1] *
                  j01, dampingConsts[2]*j01,
                  dampingConsts[0]*j02, dampingConsts[1] *
                  j02, dampingConsts[2]*j02,
                  dampingConsts[0]*j03, dampingConsts[1] *
                  j03, dampingConsts[2]*j03,
                  dampingConsts[0]*j10, dampingConsts[1] *
                  j10, dampingConsts[2]*j10,
                  dampingConsts[0]*j11, dampingConsts[1] *
                  j11, dampingConsts[2]*j11,
                  dampingConsts[0]*j12, dampingConsts[1] *
                  j12, dampingConsts[2]*j12,
                  dampingConsts[0]*j13, dampingConsts[1] *
                  j13, dampingConsts[2]*j13,
                  dampingConsts[0]*j20, dampingConsts[1] *
                  j20, dampingConsts[2]*j20,
                  dampingConsts[0]*j21, dampingConsts[1] *
                  j21, dampingConsts[2]*j21,
                  dampingConsts[0]*j22, dampingConsts[1] *
                  j22, dampingConsts[2]*j22,
                  dampingConsts[0]*j23, dampingConsts[1] *
                  j23, dampingConsts[2]*j23,
                  dampingConsts[0]*j30, dampingConsts[1] *
                  j30, dampingConsts[2]*j30,
                  dampingConsts[0]*j31, dampingConsts[1] *
                  j31, dampingConsts[2]*j31,
                  dampingConsts[0]*j32, dampingConsts[1] *
                  j32, dampingConsts[2]*j32,
                  dampingConsts[0]*j33, dampingConsts[1]*j33, dampingConsts[2]*j33]
            index = 0
            for i in range(4):
                for j in range(4):
                    jacobian[cell[i], 0, cell[j], 0] += numpy.array(jj[index])
                    jacobian[cell[i], 1, cell[j],
                             1] += numpy.array(jj[index+1])
                    jacobian[cell[i], 2, cell[j],
                             2] += numpy.array(jj[index+2])
                    index += 3

    print "Absorbing boundary for hex8 mesh"
    print "damping constants:"
    for v in dampingConsts:
        print "  %16.8e" % v
    print "values for residual:"
    for v in numpy.ravel(residual):
        print "  %16.8e," % v
    print "values for jacobian:"
    for j in numpy.ravel(jacobian):
        print "  %16.8e" % j


# ----------------------------------------------------------------------
# calcTri3()
# calcQuad4()
# calcTet4()
calcHex8()


# End of file
