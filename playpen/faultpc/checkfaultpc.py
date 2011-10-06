# Jacobian matrix for twocells/twoquad4/dislocation.cfg.

import numpy
import numpy.linalg as linalg

K = numpy.array([[4/3.0, -0.5, 1/6.0, 0,   0,0,0,0],
                 [-0.5, 4/3.0, 0, -5/6.0,  0,0,0,0],
                 [1/6.0, 0, 4/3.0, 0.5,    0,0,0,0],
                 [0, -5/6.0, 0.5, 4/3.0,   0,0,0,0],
                 [0,0,0,0,  4/3.0, 0.5, 1/6.0, 0],
                 [0,0,0,0,  0.5, 4/3.0, 0, -5/6.0],
                 [0,0,0,0, 1/6.0, 0, 4/3.0, -0.5],
                 [0,0,0,0, 0, -5/6.0, -0.5, 4/3.0]],
                dtype=numpy.float64)
Ki = linalg.inv(K)

# Gauss quadrature
L = numpy.array([[+2/3.0,0.0, +1/3.0,0.0,  -2/3.0,0.0, -1/3.0,0.0],
                 [0.0,+2/3.0, 0.0,+1/3.0,  0.0,-2/3.0, 0.0,-1/3.0],
                 [+1/3.0,0.0, +2/3.0,0.0,  -1/3.0,0.0, -2/3.0,0.0],
                 [0.0,+1/3.0, 0.0,+2/3.0,  0.0,-1/3.0, 0.0,-2/3.0]],
                dtype=numpy.float64)

# Collocated vertices and quadrature
Lc = numpy.array([[-1.0,0.0, 0.0,0.0,  +1.0,0.0, 0.0,0.0],
                  [0.0,-1.0, 0.0,0.0,  0.0,+1.0, 0.0,0.0],
                  [0.0,0.0, -1.0,0.0,  0.0,0.0, +1.0,0.0],
                  [0.0,0.0, 0.0,-1.0,  0.0,0.0, 0.0,+1.0]],
                  dtype=numpy.float64)

Z = numpy.zeros( (4,4), dtype=numpy.float64)

A = numpy.vstack( (numpy.hstack( (K, L.transpose()) ),
                   numpy.hstack( (L, Z) ) ) )
Ainv = linalg.inv(A)

Ac = numpy.vstack( (numpy.hstack( (K, Lc.transpose()) ),
                     numpy.hstack( (Lc, Z) ) ) )
Acinv = linalg.inv(Ac)

# Compute [L] [K]^(-1) [L]^T and its inverse.
LKiL = numpy.dot(numpy.dot(L, Ki), L.transpose())
LKiLi = numpy.linalg.inv(LKiL)

LKiLc = numpy.dot(numpy.dot(Lc, Ki), Lc.transpose())
LKiLic = numpy.linalg.inv(LKiLc)

# Compute diagonal approximation of LKiL and its inverse
Kid = 1.0 / K.diagonal() * numpy.identity(K.shape[0])
LKidL = numpy.dot(numpy.dot(L, Kid), L.transpose())
LKidLd = LKidL.diagonal() * numpy.identity(LKidL.shape[0])

LKidLc = numpy.dot(numpy.dot(Lc, Kid), Lc.transpose())
LKidLdc = LKidLc.diagonal() * numpy.identity(LKidLc.shape[0])

# Compute preconditioner using full matrices (no approximations)
P = A
Pi = numpy.linalg.inv(P)

# Compute condition number
evals, evecs = numpy.linalg.eig(numpy.dot(Pi, A))
print numpy.abs(evals)
print numpy.max(numpy.abs(evals))/numpy.min(numpy.abs(evals))

# Compute preconditioner using diagonal approximations for K
Pd = numpy.zeros(A.shape)
Pd[0:8,0:8] = K
Pd[0:8,8:12] = 0.0
Pd[8:12,0:8] = 0.0
Pd[8:12,8:12] = -LKidL
print LKidL

# Compute preconditioner using diagonal approximations for K
Pdc = numpy.zeros(Ac.shape)
Pdc[0:8,0:8] = K
Pdc[0:8,8:12] = 0.0
Pdc[8:12,0:8] = 0.0
Pdc[8:12,8:12] = -LKidLc
print LKidLc

# Compute condition number for diagonal approximations
Pdi = numpy.linalg.inv(Pd)
evals, evecs = numpy.linalg.eig(numpy.dot(Pdi, A))
print numpy.abs(evals)
print numpy.max(numpy.abs(evals))/numpy.min(numpy.abs(evals))

Pdci = numpy.linalg.inv(Pdc)
evals, evecs = numpy.linalg.eig(numpy.dot(Pdci, Ac))
print numpy.abs(evals)
print numpy.max(numpy.abs(evals))/numpy.min(numpy.abs(evals))
