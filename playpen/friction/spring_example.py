# Simple spring system to test friction implementation. Sliding occurs
# between DOF 2 and 3. We specify the displacement at DOF 5.
#
# |--\/\/\--*--\/\/\--*--\/\/\--| |--\/\/\--*--\/\/\--*
#      k0  u0    k1  u1    k2  u2 u3   k3  u4    k4   u5
#
# Want to solve A u(t+dt) = b(t+dt).
#
# Use displacement increment formulation, u(t+dt) = u(t) + du(t).
#
# Solve for displacement increment using A du(t) = b(t+dt) - A u(t).
#
# Residual at time t+dt is
# r(t+dt) = b(t+dt) - A (u(t) + du(t)).
#
# Constraint is the u3 - u2 = 0, unless force exceeds static friction
# (fy). If force exceeds friction, the friction is equal fy.
#
# Use Lagrange multipliers to enforce constraint.
#
# A u = b, becomes
#
# [ A C^T ] [ u ] = [ b ]
# [ C  0  ] [ l ]   [ d ]
#
# where C is the constraint matrix, l is the vector of Lagrange
# multipliers, and d is the vector of slip.

import numpy

# Spring stiffness [k0, k1, ..., kN]
k = (1.0, 1.0, 1.0, 1.0, 1.0)

# Prescribed displacement, u5
u5 = 2.5

# Static friction force
fy = 0.4

# Jacobian matrix of system
# [ K  C^T ]
# [ C  0   ]
# K is stiffness matrix
# C is the constraint matrix 
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
def reformResidual(disp, incr):
    """
    Calculate residual
    """
    calcFriction(disp, incr) # FaultCohesiveDyn::integrateResidual()
    residual = b - numpy.dot(A, disp+incr)
    residual[4] += k[4]*u5 # Dirichlet BC
    return residual


# ----------------------------------------------------------------------
def calcFriction(disp, incr):
    """
    Calculate friction.
    """
    dispTpdt = disp + incr
    lm = dispTpdt[5] # Lagrange multiplier
    if lm > fy or (lm < fy and b[5] > 0.0):
        incr[5] = fy - disp[5] # Adjust Lagrange multiplier in solution
        #slipIncr = 2*(lm-fy)
        slipIncr = (k[2] + k[3])*(lm-fy) - (k[2]*b[2] - k[3]*b[3])
        b[5] += slipIncr # Update slip estimate
        incr[2] -= 0.5*slipIncr
        incr[3] += 0.5*slipIncr


# ----------------------------------------------------------------------
def solve(incr, jacobian, residual, disp):
    """
    Solve for increment of displacment given jacobian and current disp.
    """
    iter = 0
    while numpy.sum(numpy.abs(residual)) > 1.0e-8 and iter < 40:
        print "Interation: %d" % iter
        dincr = numpy.dot(Ai, residual) # Increment to disp increment.
        incr += dincr
        residual = reformResidual(disp, incr)
        print "Disp(t):",disp
        print "Incr(t):",incr
        print "b:",b
        print "Residual:",residual
        print "dincr:",dincr
        iter += 1
    return

# ----------------------------------------------------------------------
# main
incr = 0*disp
residual = reformResidual(disp, incr)
solve(incr, A, residual, disp)
disp += incr
print "Solution:",disp
