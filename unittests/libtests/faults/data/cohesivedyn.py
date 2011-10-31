cell = "tri3d"
dim = "2d"
testCase = "open"

import numpy

from numpy import *
from numpy.linalg import inv

numpy.set_printoptions(precision=12)

# ----------------------------------------------------------------------
def printdata(data):
    """
    Print data as C array.
    """
    (nrows, ncols) = data.shape
    style = " %16.12f,"*ncols
    for row in xrange(nrows):
        print (style % tuple(data[row,:]))
    return


# ----------------------------------------------------------------------
def globalToFault(v, R):
    """
    Convert vector from global coordinate system to fault coordinate system.
    """
    (m,ndof) = v.shape

    vF = numpy.dot(C, v.reshape(m*ndof,1))
    return vF.reshape((m, ndof))


# ----------------------------------------------------------------------
def faultToGlobal(v, R):
    """
    Convert vector from fault coordinate system to global coordinate system.
    """
    (m,ndof) = v.shape

    vG = numpy.dot(C.transpose(), v.reshape(m*ndof,1))
    return vG.reshape((m, ndof))


# ----------------------------------------------------------------------
if dim == "2d":
    if cell == "tri3":
        dlagrange1 = numpy.zeros(2)
        indexL = numpy.arange(12,16)
        indexN = numpy.arange(2,6)
        indexP = numpy.arange(8,12)
        n = 16
        m = 4
        DOF = 2

        fieldT = numpy.array([[8.6, 9.6],
                              [8.8, 9.8]])
        fieldIncr = numpy.array([[1.6, 2.6],
                                 [1.8, 2.8]])
        L = numpy.array([[1.0, 0.0, 0.0, 0.0,],
                         [0.0, 1.0, 0.0, 0.0,],
                         [0.0, 0.0, 1.0, 0.0,],
                         [0.0, 0.0, 0.0, 1.0,],]);
        C = numpy.array([[0.0, -1.0, 0.0, 0.0,],
                         [-1.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, -1.0,],
                         [0.0, 0.0, -1.0, 0.0,],]);
    
        jacobianN = numpy.array(
            [[  4.0,  -1.2,  -2.2,  -2.3,],
             [  -1.2,  5.0,  -1.3,  -3.2,],
             [  -2.2,  -1.3,  4.1,  -4.3,],
             [  -2.3,  -3.2,  -4.3,  5.1,],])

        jacobianP = numpy.array(
            [[  5.0,  -1.2,  -2.2,  -2.3,],
             [  -1.2,  4.0,  -1.3,  -3.2,],
             [  -2.2,  -1.3,  5.1,  -4.3,],
             [  -2.3,  -3.2,  -4.3,  4.1,],])

        disp = numpy.array([[ 8.1, 9.1,],
                            [ 8.2, 9.2,],
                            [ 8.3, 9.3,],
                            [ 8.4, 9.4,],
                            [ 8.2, 9.2,],
                            [ 8.3, 9.3,],
                            [ 8.6, 9.6,],
                            [ 8.8, 9.8,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 9.1, 7.1,],
                                    [ 9.2, 7.2,],
                                    [ 9.3, 7.3,],
                                    [ 9.4, 7.4,],
                                    [ 9.2, 7.2,],
                                    [ 9.3, 7.3,],
                                    [ 1.6, 2.6,],
                                    [ 1.8, 2.8,],])            
        elif testCase == "open":
            dispIncr = numpy.array([[ 9.1, 7.1,],
                                    [ 9.2, 7.2,],
                                    [ 9.3, 7.3,],
                                    [ 9.4, 7.4,],
                                    [ 9.2, 7.2,],
                                    [ 9.3, 7.3,],
                                    [ -10.6, 2.6,],
                                    [ -10.8, 2.8,],])            


    elif cell == "tri3d":
        dlagrange1 = numpy.zeros(3)
        indexL = numpy.array([18, 19, 20, 21, 22, 23])
        indexN = numpy.array([2, 3, 4, 5, 8, 9])
        indexP = numpy.array([12, 13, 14, 15, 16, 17])
        n = 24
        m = 6
        DOF = 2

        fieldT = numpy.array([[-3.8, 4.8],
                              [3.0, 4.0],
                              [3.2, 4.2]])
        fieldIncr = numpy.array([[1.8, 0.8],
                                 [1.0, 0.1],
                                 [1.2, 0.2]])

        L = numpy.array([[2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])

        C = numpy.array([[+0.70710678118654757, -0.70710678118654757, 0.0, 0.0, 0.0, 0.0,],
                         [-0.70710678118654757, -0.70710678118654757, 0.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, -1.0, 0.0, 0.0,],
                         [0.0, 0.0, -1.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, 0.0, +1.0, 0.0,],
                         [0.0, 0.0, 0.0, 0.0, 0.0, -1.0,],])
    
        jacobianN = numpy.array(
            [[+6.0, -1.0, -1.1, -1.2, -1.3, -1.4],
             [-1.0, +6.1, -0.9, -0.8, -0.7, -0.6],
             [-1.1, -0.9, +6.2, -2.1,  0.0,  0.0],
             [-1.2, -0.8, -2.1, +6.3,  0.0,  0.0],
             [-1.3, -0.7,  0.0,  0.0, +6.4, -1.1],
             [-1.4, -0.6,  0.0,  0.0, -1.1, +6.5]])

        jacobianP = numpy.array(
            [[+5.0, -1.0, -1.1, -1.2, -1.3, -1.4],
             [-1.0, +5.1, -0.9, -0.8, -0.7, -0.6],
             [-1.1, -0.9, +5.2, -2.1,  0.0,  0.0],
             [-1.2, -0.8, -2.1, +5.3,  0.0,  0.0],
             [-1.3, -0.7,  0.0,  0.0, +5.4, -1.1],
             [-1.4, -0.6,  0.0,  0.0, -1.1, +5.5]])

        disp = numpy.array([[ 6.1, 8.1,],
                            [ 6.2, 8.2,],
                            [ 6.3, 8.3,],
                            [ 6.4, 8.4,],
                            [ 6.5, 8.5,],
                            [ 6.6, 8.6,],
                            [ 6.2, 8.2,],
                            [ 6.3, 8.3,],
                            [ 6.5, 8.5,],
                            [-3.8, 4.8,],
                            [ 3.0, 4.0,],
                            [ 3.2, 4.2,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 1.1, 2.1,],
                                    [ 1.2, 2.2,],
                                    [ 1.3, 2.3,],
                                    [ 1.4, 2.4,],
                                    [ 1.5, 2.5,],
                                    [ 1.6, 2.6,],
                                    [ 1.2, 2.2,],
                                    [ 1.3, 2.3,],
                                    [ 1.5, 2.5,],
                                    [ 1.8, 0.8,],
                                    [ 1.0, 0.1,],
                                    [ 1.2, 0.2,],])            
        elif testCase == "open":
            dispIncr = numpy.array([[ 1.1, 2.1,],
                                    [ 1.2, 2.2,],
                                    [ 1.3, 2.3,],
                                    [ 1.4, 2.4,],
                                    [ 1.5, 2.5,],
                                    [ 1.6, 2.6,],
                                    [ 1.2, 2.2,],
                                    [ 1.3, 2.3,],
                                    [ 1.5, 2.5,],
                                    [-10.8, 0.8,],
                                    [-10.0, 0.1,],
                                    [ 1.2, -10.2,],])            


    elif cell == "quad4":
        dlagrange1 = numpy.zeros(2)
        indexL = numpy.arange(16,20)
        indexN = numpy.arange(4,8)
        indexP = numpy.arange(12,16)
        n = 20
        m = 4
        DOF = 2

        fieldT = numpy.array([[8.8, 9.8],
                              [8.0, 9.0]])
        fieldIncr = numpy.array([[1.8, 2.8],
                                 [1.0, 2.0]])
        L = numpy.array([[1.0, 0.0, 0.0, 0.0,],
                         [0.0, 1.0, 0.0, 0.0,],
                         [0.0, 0.0, 1.0, 0.0,],
                         [0.0, 0.0, 0.0, 1.0,],]);
        C = numpy.array([[0.0, -1.0, 0.0, 0.0,],
                         [-1.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, -1.0,],
                         [0.0, 0.0, -1.0, 0.0,],]);
    
        jacobianN = numpy.array(
            [[  4.0,  -1.2,  -2.2,  -2.3,],
             [  -1.2,  5.0,  -1.3,  -3.2,],
             [  -2.2,  -1.3,  4.1,  -4.3,],
             [  -2.3,  -3.2,  -4.3,  5.1,],])

        jacobianP = numpy.array(
            [[  5.0,  -1.2,  -2.2,  -2.3,],
             [  -1.2,  4.0,  -1.3,  -3.2,],
             [  -2.2,  -1.3,  5.1,  -4.3,],
             [  -2.3,  -3.2,  -4.3,  4.1,],])

        disp = numpy.array([[ 8.1, 9.1,],
                            [ 8.2, 9.2,],
                            [ 8.3, 9.3,],
                            [ 8.4, 9.4,],
                            [ 8.5, 9.5,],
                            [ 8.6, 9.6,],
                            [ 8.3, 9.3,],
                            [ 8.4, 9.4,],
                            [ 8.8, 9.8,],
                            [ 8.0, 9.0,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 1.1, 2.1,],
                                    [ 1.2, 2.2,],
                                    [ 1.3, 2.3,],
                                    [ 1.4, 2.4,],
                                    [ 1.5, 2.5,],
                                    [ 1.6, 2.6,],
                                    [ 1.3, 2.3,],
                                    [ 1.4, 2.4,],
                                    [ 1.8, 2.8,],
                                    [ 1.0, 2.0,],])            
        elif testCase == "open":
            dispIncr = numpy.array([[ 1.1, 2.1,],
                                    [ 1.2, 2.2,],
                                    [ 1.3, 2.3,],
                                    [ 1.4, 2.4,],
                                    [ 1.5, 2.5,],
                                    [ 1.6, 2.6,],
                                    [ 1.3, 2.3,],
                                    [ 1.4, 2.4,],
                                    [ -10.8, 2.8,],
                                    [ -10.0, 2.0,],])            


    # ------------------------------------------------------------------
    fieldTpdt = fieldT + fieldIncr

    fieldTpdt = globalToFault(fieldTpdt, C)

    tractionShear = abs(fieldTpdt[:,0])
    tractionNormal = fieldTpdt[:,1]

    print "tractionShear",tractionShear
    print "tractionNormal",tractionNormal

    friction = -0.6 * tractionNormal;

    print "friction",friction

    dlagrange0 = (friction - tractionShear) * fieldTpdt[:,0] / tractionShear
                           
    print "dlagrange0",dlagrange0

    if testCase == "slip": 
        dLagrange = numpy.vstack((dlagrange0, dlagrange1))
        dLagrange = numpy.transpose(dLagrange)
        dLagrange = faultToGlobal(dLagrange, C).reshape(m)
    elif testCase == "open":
        dLagrange = numpy.reshape(disp+dispIncr, n)
        dLagrange = -dLagrange[indexL]

    print "dLagrange \n", dLagrange

    RHS = numpy.dot(numpy.transpose(L),dLagrange)
    print "RHS",RHS
    duN = numpy.dot(inv(jacobianN),RHS)
    duP = -numpy.dot(inv(jacobianP),RHS)
    
    dispRel = duP - duN

    dispTpdt = disp + dispIncr
    dispTpdt = numpy.reshape(dispTpdt, n)

    slipVertex = dispRel + dispTpdt[indexP]-dispTpdt[indexN]
    slipVertex = numpy.reshape(slipVertex, (m/DOF,DOF))
    slipVertex = globalToFault(slipVertex, C)
    if testCase == "slip":
        slipVertex[:,1] = 0
    mask = slipVertex[:,1] < 0.0
    slipVertex[mask,1] = 0
    print "slip",slipVertex
    slipVertex = faultToGlobal(slipVertex, C)
    slipVertex = numpy.reshape(slipVertex, m)

    print "duN \n", duN
    print "duP \n", duP

    dispIncrE = dispIncr
    dispIncrE = numpy.reshape(dispIncrE, n)
    dispIncrE[indexL] = dispIncrE[indexL] + dLagrange
    dispIncrE[indexN] = dispIncrE[indexN] - 0.5*slipVertex
    dispIncrE[indexP] = dispIncrE[indexP] + 0.5*slipVertex
    dispIncrE = numpy.reshape(dispIncrE, (n/DOF,DOF))

    slipVertex = numpy.reshape(slipVertex, (m/DOF,DOF))
    slipVertex = globalToFault(slipVertex, C)

    print "dispIncrE\n", printdata(dispIncrE)
    print "slipVertexE\n", printdata(slipVertex)


# ----------------------------------------------------------------------
elif dim == "3d":
    if cell == "tet4":

        dlagrange2 = numpy.zeros(3)
        indexL = numpy.arange(24,33)
        indexN = numpy.arange(3,12)
        indexP = numpy.arange(15,24)
        n = 33
        m = 9
        DOF = 3

        fieldT = numpy.array([[7.7, 8.7, 9.7],
                              [7.9, 8.9, 9.9],
                              [7.1, 8.1, 9.1]])
        fieldIncr = numpy.array([[9.7, 2.7, 3.7],
                                 [9.9, 2.9, 3.9],
                                 [9.1, 2.1, 3.1]])
        
        L = numpy.array([[1.0/3.0,0,0, 0.0,0,0, 0.0,0,0,],
                         [0,1.0/3.0,0, 0,0.0,0, 0,0.0,0,],
                         [0,0,1.0/3.0, 0,0,0.0, 0,0,0.0,],
                         [0.0,0,0, 1.0/3.0,0,0, 0.0,0,0,],
                         [0,0.0,0, 0,1.0/3.0,0, 0,0.0,0,],
                         [0,0,0.0, 0,0,1.0/3.0, 0,0,0.0,],
                         [0.0,0,0, 0.0,0,0, 1.0/3.0,0,0,],
                         [0,0.0,0, 0,0.0,0, 0,1.0/3.0,0,],
                         [0,0,0.0, 0,0,0.0, 0,0,1.0/3.0,]])

        Cv = numpy.array([[ 0, -1, 0,],
                          [ 0, 0, +1,],
                          [ -1, 0, 0,],])
        Zv = numpy.zeros([3,3])
        C = numpy.vstack( (numpy.hstack((Cv, Zv, Zv)),
                           numpy.hstack((Zv, Cv, Zv)),
                           numpy.hstack((Zv, Zv, Cv)) ) )

        jacobianN = numpy.array(
            [[ 4.0, -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8],
             [-1.1,  4.1, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9], 
             [-1.2, -2.3,  4.2, -1.0, -1.1, -1.2, -1.3, -1.4, -1.5],
             [-1.3, -2.4, -1.0,  4.3, -0.2, -0.3, -0.4, -0.5, -0.6],
             [-1.4, -2.5, -1.1, -0.2,  4.4, -0.9, -0.8, -0.7, -0.5],
             [-1.5, -2.6, -1.2, -0.3, -0.9,  4.5, -1.1, -1.2, -1.3],
             [-1.6, -2.7, -1.3, -0.4, -0.8, -1.1,  4.6, -1.8, -1.5],
             [-1.7, -2.8, -1.4, -0.5, -0.7, -1.2, -1.8,  4.7, -1.1],
             [-1.8, -2.9, -1.5, -0.6, -0.5, -1.3, -1.5, -1.1,  4.8]])

        jacobianP = numpy.array(
            [[ 5.0, -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8],
             [-1.1,  5.1, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9], 
             [-1.2, -2.3,  5.2, -1.0, -1.1, -1.2, -1.3, -1.4, -1.5],
             [-1.3, -2.4, -1.0,  5.3, -0.2, -0.3, -0.4, -0.5, -0.6],
             [-1.4, -2.5, -1.1, -0.2,  5.4, -0.9, -0.8, -0.7, -0.5],
             [-1.5, -2.6, -1.2, -0.3, -0.9,  5.5, -1.1, -1.2, -1.3],
             [-1.6, -2.7, -1.3, -0.4, -0.8, -1.1,  5.6, -1.8, -1.5],
             [-1.7, -2.8, -1.4, -0.5, -0.7, -1.2, -1.8,  5.7, -1.1],
             [-1.8, -2.9, -1.5, -0.6, -0.5, -1.3, -1.5, -1.1,  5.8]])

        disp = numpy.array([[ 7.1, 8.1, 9.1,],
                            [ 7.2, 8.2, 9.2,],
                            [ 7.3, 8.3, 9.3,],
                            [ 7.4, 8.4, 9.4,],
                            [ 7.5, 8.5, 9.5,],
                            [ 7.2, 8.2, 9.2,],
                            [ 7.3, 8.3, 9.3,],
                            [ 7.4, 8.4, 9.4,],
                            [ 7.7, 8.7, 9.7,],
                            [ 7.9, 8.9, 9.9,],
                            [ 7.1, 8.1, 9.1,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 1.1, 2.1, 3.1,],
                                    [ 1.2, 2.2, 3.2,],
                                    [ 1.3, 2.3, 3.3,],
                                    [ 1.4, 2.4, 3.4,],
                                    [ 1.5, 2.5, 3.5,],
                                    [ 1.2, 2.2, 3.2,],
                                    [ 1.3, 2.3, 3.3,],
                                    [ 1.4, 2.4, 3.4,],
                                    [ 9.7, 2.7, 3.7,],
                                    [ 9.9, 2.9, 3.9,],
                                    [ 9.1, 2.1, 3.1,],])            
        elif testCase == "open":
            dispIncr = numpy.array([[ 1.1, 2.1, 3.1,],
                                    [ 1.2, 2.2, 3.2,],
                                    [ 1.3, 2.3, 3.3,],
                                    [ 1.4, 2.4, 3.4,],
                                    [ 1.5, 2.5, 3.5,],
                                    [ 1.2, 2.2, 3.2,],
                                    [ 1.3, 2.3, 3.3,],
                                    [ 1.4, 2.4, 3.4,],
                                    [-20.7,  2.7, 3.7,],
                                    [-20.9,  2.9, 3.9,],
                                    [-20.1,  2.1, 3.1,],])            


    elif cell == "hex8":
        dlagrange2 = numpy.zeros(4)
        indexL = numpy.arange(48,60)
        indexN = numpy.arange(12,24)
        indexP = numpy.arange(36,48)
        n = 60
        m = 12
        DOF = 3

        a0 = 1.0
        a1 = 0.0
        a2 = 0.0
        L = numpy.array([[a0, 0, 0, a1, 0, 0, a1, 0, 0, a2, 0, 0],
                         [0, a0, 0, 0, a1, 0, 0, a1, 0, 0, a2, 0],
                         [0, 0, a0, 0, 0, a1, 0, 0, a1, 0, 0, a2],
                         [a1, 0, 0, a0, 0, 0, a2, 0, 0, a1, 0, 0],
                         [0, a1, 0, 0, a0, 0, 0, a2, 0, 0, a1, 0],
                         [0, 0, a1, 0, 0, a0, 0, 0, a2, 0, 0, a1],
                         [a1, 0, 0, a2, 0, 0, a0, 0, 0, a1, 0, 0],
                         [0, a1, 0, 0, a2, 0, 0, a0, 0, 0, a1, 0],
                         [0, 0, a1, 0, 0, a2, 0, 0, a0, 0, 0, a1],
                         [a2, 0, 0, a1, 0, 0, a1, 0, 0, a0, 0, 0],
                         [0, a2, 0, 0, a1, 0, 0, a1, 0, 0, a0, 0],
                         [0, 0, a2, 0, 0, a1, 0, 0, a1, 0, 0, a0]])

        fieldT = numpy.array([[4.4, 2.4, 3.4],
                              [4.6, 2.6, 3.6],
                              [4.8, 2.8, 3.8],
                              [4.0, 2.0, 3.0]])
        fieldIncr = numpy.array([[1.4, 2.4, 0.4],
                                 [1.6, 2.6, 0.6],
                                 [1.8, 2.8, 0.8],
                                 [1.0, 2.0, 0.2]])

        Cv = numpy.array([[ 0, -1, 0,],
                          [ 0, 0, +1,],
                          [ -1, 0, 0,],])
        Zv = numpy.zeros([3,3])
        C = numpy.vstack( (numpy.hstack((Cv, Zv, Zv, Zv)),
                           numpy.hstack((Zv, Cv, Zv, Zv)),
                           numpy.hstack((Zv, Zv, Cv, Zv)),
                           numpy.hstack((Zv, Zv, Zv, Cv)) ) )

        jacobianN = numpy.array(
            [[+6.0, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0, -0.8, -0.7, -0.6, -0.5, -0.4,],
             [-0.5, +6.1, -1.0, -1.1, -1.2, -1.3, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9,],
             [-0.6, -1.0, +6.2, -0.5, -0.6, -0.7, -0.8, -0.9, -0.8, -0.7, -0.6, -0.5,],
             [-0.7, -1.1, -0.5, +6.3, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,],
             [-0.8, -1.2, -0.6, -0.8, +6.4, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9,],
             [-0.9, -1.3, -0.7, -0.7, -0.3, +6.5, -0.3, -0.8, -0.7, -0.6, -0.9, -0.7,],
             [-1.0, -1.4, -0.8, -0.6, -0.4, -0.3, +6.6, -1.1, -0.8, -0.7, -0.6, -0.5,],
             [-0.8, -1.3, -0.9, -0.5, -0.5, -0.8, -1.1, +6.7, -0.8, -0.9, -1.0, -1.1,],
             [-0.7, -1.2, -0.8, -0.4, -0.6, -0.7, -0.8, -0.8, +6.8, -1.0, -1.1, -1.2,],
             [-0.6, -1.1, -0.7, -0.3, -0.7, -0.6, -0.7, -0.9, -1.0, +6.9, -0.5, -0.4,],
             [-0.5, -1.0, -0.6, -0.2, -0.8, -0.9, -0.6, -1.0, -1.1, -0.5, +6.0, -1.2,],
             [-0.4, -0.9, -0.5, -0.1, -0.9, -0.7, -0.5, -1.1, -1.2, -0.4, -1.2, +6.1,],])


        jacobianP = numpy.array(
            [[+7.0, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0, -0.8, -0.7, -0.6, -0.5, -0.4,],
             [-0.5, +7.1, -1.0, -1.1, -1.2, -1.3, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9,],
             [-0.6, -1.0, +7.2, -0.5, -0.6, -0.7, -0.8, -0.9, -0.8, -0.7, -0.6, -0.5,],
             [-0.7, -1.1, -0.5, +7.3, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,],
             [-0.8, -1.2, -0.6, -0.8, +7.4, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9,],
             [-0.9, -1.3, -0.7, -0.7, -0.3, +7.5, -0.3, -0.8, -0.7, -0.6, -0.9, -0.7,],
             [-1.0, -1.4, -0.8, -0.6, -0.4, -0.3, +7.6, -1.1, -0.8, -0.7, -0.6, -0.5,],
             [-0.8, -1.3, -0.9, -0.5, -0.5, -0.8, -1.1, +7.7, -0.8, -0.9, -1.0, -1.1,],
             [-0.7, -1.2, -0.8, -0.4, -0.6, -0.7, -0.8, -0.8, +7.8, -1.0, -1.1, -1.2,],
             [-0.6, -1.1, -0.7, -0.3, -0.7, -0.6, -0.7, -0.9, -1.0, +7.9, -0.5, -0.4,],
             [-0.5, -1.0, -0.6, -0.2, -0.8, -0.9, -0.6, -1.0, -1.1, -0.5, +7.0, -1.2,],
             [-0.4, -0.9, -0.5, -0.1, -0.9, -0.7, -0.5, -1.1, -1.2, -0.4, -1.2, +7.1,],])

        disp = numpy.array([[ 4.1, 2.1, 3.1,],
                            [ 4.2, 2.2, 3.2,],
                            [ 4.3, 2.3, 3.3,],
                            [ 4.4, 2.4, 3.4,],
                            [ 4.5, 2.5, 3.5,],
                            [ 4.6, 2.6, 3.6,],
                            [ 4.7, 2.7, 3.7,],
                            [ 4.8, 2.8, 3.8,],
                            [ 4.9, 2.9, 3.9,],
                            [ 4.0, 2.0, 3.0,],
                            [ 4.1, 2.1, 3.1,],
                            [ 4.2, 2.2, 3.2,],
                            [ 4.5, 2.5, 3.5,],
                            [ 4.6, 2.6, 3.6,],
                            [ 4.7, 2.7, 3.7,],
                            [ 4.8, 2.8, 3.8,],
                            [ 4.4, 2.4, 3.4,],
                            [ 4.6, 2.6, 3.6,],
                            [ 4.8, 2.8, 3.8,],
                            [ 4.0, 2.0, 3.0,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 1.1, 2.1, 0.1,],
                                    [ 1.2, 2.2, 0.2,],
                                    [ 1.3, 2.3, 0.3,],
                                    [ 1.4, 2.4, 0.4,],
                                    [ 1.5, 2.5, 0.5,],
                                    [ 1.6, 2.6, 0.6,],
                                    [ 1.7, 2.7, 0.7,],
                                    [ 1.8, 2.8, 0.8,],
                                    [ 1.9, 2.9, 0.9,],
                                    [ 1.0, 2.0, 0.0,],
                                    [ 1.1, 2.1, 0.1,],
                                    [ 1.2, 2.2, 0.2,],
                                    [ 1.5, 2.5, 0.5,],
                                    [ 1.6, 2.6, 0.6,],
                                    [ 1.7, 2.7, 0.7,],
                                    [ 1.8, 2.8, 0.8,],
                                    [ 1.4, 2.4, 0.4,],
                                    [ 1.6, 2.6, 0.6,],
                                    [ 1.8, 2.8, 0.8,],
                                    [ 1.0, 2.0, 0.2,],])          
        elif testCase == "open":
            dispIncr = numpy.array([[ 1.1, 2.1, 0.1,],
                                    [ 1.2, 2.2, 0.2,],
                                    [ 1.3, 2.3, 0.3,],
                                    [ 1.4, 2.4, 0.4,],
                                    [ 1.5, 2.5, 0.5,],
                                    [ 1.6, 2.6, 0.6,],
                                    [ 1.7, 2.7, 0.7,],
                                    [ 1.8, 2.8, 0.8,],
                                    [ 1.9, 2.9, 0.9,],
                                    [ 1.0, 2.0, 0.0,],
                                    [ 1.1, 2.1, 0.1,],
                                    [ 1.2, 2.2, 0.2,],
                                    [ 1.5, 2.5, 0.5,],
                                    [ 1.6, 2.6, 0.6,],
                                    [ 1.7, 2.7, 0.7,],
                                    [ 1.8, 2.8, 0.8,],
                                    [-10.4, 2.4, 0.4,],
                                    [-10.6, 2.6, 0.6,],
                                    [-10.8, 2.8, 0.8,],
                                    [-10.0, 2.0, 0.2,],])          

    # ------------------------------------------------------------------
    fieldTpdt = fieldT + fieldIncr

    fieldTpdt = globalToFault(fieldTpdt, C)

    tractionShear = (fieldTpdt[:,0]**2 + fieldTpdt[:,1]**2)**0.5
    tractionNormal = fieldTpdt[:,2]

    print "tractionShear",tractionShear
    print "tractionNormal",tractionNormal

    friction = -0.6 * tractionNormal;

    print "friction",friction

    dlagrange0 = (friction - tractionShear) * fieldTpdt[:,0] / tractionShear
    dlagrange1 = (friction - tractionShear) * fieldTpdt[:,1] / tractionShear
                           
    print "dlagrange0",dlagrange0
    print "dlagrange1",dlagrange1

    if testCase == "slip": 
        dLagrange = numpy.vstack((dlagrange0, dlagrange1, dlagrange2))
        dLagrange = numpy.transpose(dLagrange)
        dLagrange = faultToGlobal(dLagrange, C).reshape(m)
    elif testCase == "open":
        dLagrange = numpy.reshape(disp+dispIncr, n)
        dLagrange = -dLagrange[indexL]

    print "dLagrange \n", dLagrange

    RHS = numpy.dot(numpy.transpose(L),dLagrange)
    print "RHS",RHS
    duN = numpy.dot(inv(jacobianN),RHS)
    duP = -numpy.dot(inv(jacobianP),RHS)
    
    dispRel = duP - duN

    dispTpdt = disp + dispIncr
    dispTpdt = numpy.reshape(dispTpdt, n)

    slipVertex = dispRel + dispTpdt[indexP]-dispTpdt[indexN]
    slipVertex = numpy.reshape(slipVertex, (m/DOF,DOF))
    slipVertex = globalToFault(slipVertex, C)
    if testCase == "slip":
        slipVertex[:,2] = 0
    mask = slipVertex[:,2] < 0.0
    slipVertex[mask,2] = 0
    slipVertex = faultToGlobal(slipVertex, C)
    slipVertex = numpy.reshape(slipVertex, m)

    print "duN \n", duN
    print "duP \n", duP

    dispIncrE = dispIncr
    dispIncrE = numpy.reshape(dispIncrE, n)

    dispIncrE[indexL] = dispIncrE[indexL] + dLagrange
    dispIncrE[indexN] = dispIncrE[indexN] - 0.5*slipVertex
    dispIncrE[indexP] = dispIncrE[indexP] + 0.5*slipVertex
    dispIncrE = numpy.reshape(dispIncrE, (n/DOF,DOF))

    slipVertex = numpy.reshape(slipVertex, (m/DOF,DOF))
    slipVertex = globalToFault(slipVertex, C)

    print "dispIncrE\n", printdata(dispIncrE)
    print "slipVertexE\n", printdata(slipVertex)
