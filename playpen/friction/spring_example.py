#
#
#
#

import numpy

# Spring stiffness
k = (1.0, 1.0, 1.0, 1.0, 1.0)

# Prescribed displacement
u0 = 2.5

# Static friction force
fy = 0.4

# Jacobian matrix of system
# [ K  C^T ]
# [ C  0   ]
A = numpy.array([[k[0]+k[1], -k[1], 0, 0, 0, 0],
                 [-k[1], k[1]+k[2], -k[2], 0, 0, 0],
                 [0, -k[2], k[2], 0, 0, -1],
                 [0, 0, 0, k[3], -k[3], 1],
                 [0, 0, 0, -k[3], k[3]+k[4], 0],
                 [0, 0, -1, 1, 0, 0]],
                dtype=numpy.float64)
disp = numpy.zeros( (6,1), dtype=numpy.float64)
b = numpy.zeros( (6,1), dtype=numpy.float64)
Ai = numpy.linalg.inv(A)

# ----------------------------------------------------------------------
def reformResidual(disp):
    """
    Calculate residual
    """
    residual = b - numpy.dot(A, disp)
    residual[4] += k[4]*u0
    return residual


# ----------------------------------------------------------------------
def calcFriction(disp, incr):
    """
    Calculate friction.
    """
    global slipPrev
    global s
    global sset
    dispTpdt = disp+incr
    lm = dispTpdt[5]
    if lm > fy or (lm < fy and b[5] > 0.0):
        incr[5] = fy - disp[5] # Adjust Lagrange multiplier in solution
        slip = 2*(lm-fy)
        b[5] += slip # Update slip estimate
        incr[2] -= 0.5*slip
        incr[3] += 0.5*slip


# ----------------------------------------------------------------------
def solve(incr, jacobian, residual, disp):
    """
    Solve for increment of displacment given jacobian and current disp.
    """
    i = 0
    while numpy.sum(numpy.abs(residual)) > 1.0e-6 and i < 40:
        print "Interation: %d" % i
        ii = numpy.dot(Ai, residual)
        incr += ii
        calcFriction(disp, incr)
        residual = reformResidual(disp+incr)
        print "Disp(t):",disp
        print "Incr(t):",incr
        print "b:",b
        print "Residual:",residual
        print "ii:",ii
        i += 1
    return

# ----------------------------------------------------------------------
# main
incr = 0*disp
residual = reformResidual(disp+incr)
solve(incr, A, residual, disp)
disp += incr
print "Solution:",disp
print A
