#!/usr/bin/env python

import numpy

ka = 1.0
kb = 1.0

def friction(slip):
    fs = 1.0
    fd = 0.2
    d0 = 1.0e-6

    f = fs;
    if slip > d0:
        f = fd
    elif slip > 0.0:
        f = fs - (fs-fd) * abs(slip) / d0;
    return f


def computeResidual(u):
    l = u[4]
    slip = u[2] - u[1]
    fc = friction(slip) - l
    d = 0.0
    u4 = 4.1

    r = numpy.zeros(u.shape)
    r[0] = -2*ka*u[0] + ka*u[1];
    r[1] = +ka*u[0] - ka*u[1] + fc;
    r[2] = -kb*u[2] + kb*u[3] - fc;
    r[3] = +kb*u[2] - 2*kb*u[3] + kb*u4;
    r[4] = l * (u[2] - u[1] - d);
    return r

    
def jacobianInv(u):
    l = u[4]
    if l > 0.0:
        df = 0.0
    else:
        df = 0.0

    slip = u[2]-u[1]
    if l > 0.0 or abs(slip) > 0:
        jl = slip
    else:
        jl = 1.0
        
    J = numpy.array([
        [-2.0,  1.0,  0.0,  0.0,  0.0],
        [ 1.0, -1.0-df, +df,  0.0, -1.0],
        [ 0.0, +df, -1.0-df,  1.0,  1.0],
        [ 0.0,  0.0,  1.0, -2.0,  0.0],
        [ 0.0,  -l,  +l,  0.0,  jl]
        ])
    print(J)
    Jinv = numpy.linalg.inv(J)
    return Jinv

    
def computeUpdate(u,Jinv):
    r = computeResidual(u)
    du = -numpy.dot(Jinv, r)
    print("Residual",r)
    print("Update",du)
    return u + du
    

u0 = numpy.array([1.0, 2.0, 2.0, 3.0, 0.0]).T
Jinv = jacobianInv(u0)

u1 = numpy.array([1.0, 2.0, 2.1, 3.1, 0.0]).T
u2 = computeUpdate(u1, Jinv)
print("New Soln",u2)
print("Residual for new soln",computeResidual(u2))
