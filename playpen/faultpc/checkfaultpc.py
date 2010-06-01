# Jacobian matrix for twocells/twoquad4/dislocation.cfg.

import numpy
import numpy.linalg as linalg

A = numpy.array([[4/3.0, -0.5, 1/6.0, 0,   0,0,0,0],
                 [-0.5, 4/3.0, 0, -5/6.0,  0,0,0,0],
                 [1/6.0, 0, 4/3.0, 0.5,    0,0,0,0],
                 [0, -5/6.0, 0.5, 4/3.0,   0,0,0,0],
                 [0,0,0,0,  4/3.0, 0.5, 1/6.0, 0],
                 [0,0,0,0,  0.5, 4/3.0, 0, -5/6.0],
                 [0,0,0,0, 1/6.0, 0, 4/3.0, -0.5],
                 [0,0,0,0, 0, -5/6.0, -0.5, 4/3.0]],
                dtype=numpy.float64)
Ai = linalg.inv(A)

C = numpy.array([[0,-1,  0,0,  0,1, 0,0],
                 [-1,0,  0,0,  1,0, 0,0],
                 [0,0,  0,-1,  0,0, 0,1],
                 [0,0, -1,0,   0,0, 1,0]],
                dtype=numpy.float64)
Z = numpy.zeros( (4,4), dtype=numpy.float64)

J = numpy.vstack( (numpy.hstack( (A, C.transpose()) ),
                   numpy.hstack( (C, Z) ) ) )
Jinv = linalg.inv(J)

# Compute [C] [A]^(-1) [C]^T and its inverse.
CAC = numpy.dot(numpy.dot(C, Ai), C.transpose())
CACi = numpy.linalg.inv(CAC)

# Compute diagonal approximation of CAC and its inverse
Aid = 1.0 / A.diagonal() * numpy.identity(A.shape[0])
CACd = numpy.dot(numpy.dot(C, Aid), C.transpose())
CACdi = 1.0 / CACd.diagonal() * numpy.identity(CACd.shape[0])

# Compute preconditioner using full matrices (no approximations)
P = J
Pi = numpy.linalg.inv(P)

# Compute condition number
evals, evecs = numpy.linalg.eig(numpy.dot(J, Pi))
print numpy.abs(evals)
print numpy.max(numpy.abs(evals))/numpy.min(numpy.abs(evals))

# Compute preconditioner using diagonal approximations (but full A)
Pd = numpy.zeros(J.shape)
Pd[0:8,0:8] = A
Pd[0:8,8:12] = 0.0
Pd[8:12,0:8] = 0.0
Pd[8:12,8:12] = -CACd

Pdi = numpy.linalg.inv(Pd)

# Compute condition number for diagonal approximations
evals, evecs = numpy.linalg.eig(numpy.dot(Pdi, J))
print numpy.abs(evals)
print numpy.max(numpy.abs(evals))/numpy.min(numpy.abs(evals))
