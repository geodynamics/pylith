cell = "tri3"
dim = "2d"
testCase = "slip"

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
        fieldIncr = numpy.array([[9.6, -10.6],
                                 [9.8, -10.8]])
        area = numpy.array([1.0, 1.0])
        C = numpy.array([[0.0, -1.0, 0.0, 0.0,],
                         [-1.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, -1.0,],
                         [0.0, 0.0, -1.0, 0.0,],]);
    
        jacobianN = numpy.array(
            [[  1.0,  5.2,  4.2,  5.3,],
             [  6.2,  1.0,  6.3,  7.2,],
             [  8.2,  9.2,  1.0,  9.3,],
             [ 10.2, 11.2, 10.3,  1.0,],])

        jacobianP = numpy.array(
            [[  1.0, 17.5, 16.5, 17.6,],
             [ 18.5,  1.0, 18.6, 19.5,],
             [ 20.5, 21.5,  1.0, 21.6,],
             [ 22.5, 23.5, 22.6,  1.0,],])

        disp = numpy.array([[ 8.1, 9.1,],
                            [ 8.2, 9.2,],
                            [ 8.3, 9.3,],
                            [ 8.4, 9.4,],
                            [ 8.5, 9.5,],
                            [ 8.7, 9.7,],
                            [ 8.6, 9.6,],
                            [ 8.8, 9.8,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 9.1, 10.1,],
                                    [ 9.2, 10.2,],
                                    [ 9.3, 10.3,],
                                    [ 9.4, 10.4,],
                                    [ 9.5, 10.5,],
                                    [ 9.7, 10.7,],
                                    [ 9.6, -10.6,],
                                    [ 9.8, -10.8,],])            
        elif testCase == "open":
            dispIncr = numpy.array([[ 9.1, 10.1,],
                                    [ 9.2, 10.2,],
                                    [ 9.3, 10.3,],
                                    [ 9.4, 10.4,],
                                    [ 9.5, 10.5,],
                                    [ 9.7, 10.7,],
                                    [ 9.6, 10.6,],
                                    [ 9.8, 10.8,],])


    elif cell == "tri3d":
        dlagrange1 = numpy.zeros(3)
        indexL = numpy.array([18, 19, 20, 21, 22, 23])
        indexN = numpy.array([2, 3, 4, 5, 8, 9])
        indexP = numpy.array([12, 13, 14, 15, 16, 17])
        n = 24
        m = 6
        DOF = 2

        fieldT = numpy.array([[6.8, 8.8],
                              [6.0, 8.0],
                              [7.2, 9.2]])
        fieldIncr = numpy.array([[9.8, -10.8],
                                 [9.0, -10.0],
                                 [9.2, -10.2]])
        area = numpy.array([2.0, 1.0, 1.0])
        C = numpy.array([[+0.70710678118654757, -0.70710678118654757, 0.0, 0.0, 0.0, 0.0,],
                         [-0.70710678118654757, -0.70710678118654757, 0.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, -1.0, 0.0, 0.0,],
                         [0.0, 0.0, -1.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, 0.0, +1.0, 0.0,],
                         [0.0, 0.0, 0.0, 0.0, 0.0, -1.0,],])
    
        jacobianN = numpy.array(
            [[6.0, 7.3, 6.9, 7.4, 6.7, 7.6],
             [3.1, 5.2, 3.2, 5.2, 3.4, 5.3],
             [3.9, 2.7, 3.8, 2.8, 0.0, 0.0],
             [4.1, 6.6, 4.2, 6.2, 0.0, 0.0],
             [7.9, 8.1, 0.0, 0.0, 7.4, 8.5],
             [6.4, 3.4, 0.0, 0.0, 6.1, 3.8]])

        jacobianP = numpy.array(
            [[1.6, 4.1, 1.7, 4.2, 1.8, 4.3],
             [4.6, 4.8, 4.7, 4.4, 4.8, 4.2],
             [6.9, 5.3, 7.0, 5.9, 0.0, 0.0],
             [7.2, 6.6, 7.3, 6.5, 0.0, 0.0],
             [8.4, 7.8, 0.0, 0.0, 8.3, 7.1],
             [6.3, 8.6, 0.0, 0.0, 4.7, 8.7]])

        disp = numpy.array([[ 6.1, 8.1,],
                            [ 6.2, 8.2,],
                            [ 6.3, 8.3,],
                            [ 6.4, 8.4,],
                            [ 6.5, 8.5,],
                            [ 6.6, 8.6,],
                            [ 6.7, 8.7,],
                            [ 6.9, 8.9,],
                            [ 7.1, 9.1,],
                            [ 6.8, 8.8,],
                            [ 6.0, 8.0,],
                            [ 7.2, 9.2,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 9.1, 10.1,],
                                    [ 9.2, 10.2,],
                                    [ 9.3, 10.3,],
                                    [ 9.4, 10.4,],
                                    [ 9.5, 10.5,],
                                    [ 9.6, 10.6,],
                                    [ 9.7, 10.7,],
                                    [ 9.9, 10.9,],
                                    [ 9.1, 10.1,],
                                    [ 9.8, -10.8,],
                                    [ 9.0, -10.0,],
                                    [ 9.2, -10.2,],])            
        elif testCase == "open":
            dispIncr = numpy.array([[ 9.1, 10.1,],
                                    [ 9.2, 10.2,],
                                    [ 9.3, 10.3,],
                                    [ 9.4, 10.4,],
                                    [ 9.5, 10.5,],
                                    [ 9.6, 10.6,],
                                    [ 9.7, 10.7,],
                                    [ 9.9, 10.9,],
                                    [ 9.1, 10.1,],
                                    [ 9.8, 10.8,],
                                    [ 9.0, 10.0,],
                                    [ 9.2, 10.2,],])            


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
        fieldIncr = numpy.array([[-9.8, -10.8],
                                 [-9.0, -10.0]])
        area = numpy.array([1.0, 1.0])
        C = numpy.array([[0.0, -1.0, 0.0, 0.0,],
                         [-1.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, -1.0,],
                         [0.0, 0.0, -1.0, 0.0,],]);
    
        jacobianN = numpy.array(
            [[  1.0,  8.1,  8.2,  8.3,],
             [  9.9,  1.0, 10.0, 10.1,],
             [ 11.7, 11.8,  1.0, 11.9,],
             [ 13.5, 13.6, 13.7,  1.0,],])

        jacobianP = numpy.array(
            [[  1.0, 23.7, 23.8, 23.9,],
             [ 25.5,  1.0, 25.6, 25.7,],
             [ 27.3, 27.4,  1.0, 27.5,],
             [ 29.1, 29.2, 29.3,  1.0,],])

        disp = numpy.array([[ 8.1, 9.1,],
                            [ 8.2, 9.2,],
                            [ 8.3, 9.3,],
                            [ 8.4, 9.4,],
                            [ 8.5, 9.5,],
                            [ 8.6, 9.6,],
                            [ 8.7, 9.7,],
                            [ 8.9, 9.9,],
                            [ 8.8, 9.8,],
                            [ 8.0, 9.0,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 9.1, 10.1,],
                                    [ 9.2, 10.2,],
                                    [ 9.3, 10.3,],
                                    [ 9.4, 10.4,],
                                    [ 9.5, 10.5,],
                                    [ 9.6, 10.6,],
                                    [ 9.7, 10.7,],
                                    [ 9.9, 10.9,],
                                    [ -9.8, -10.8,],
                                    [ -9.0, -10.0,],])
          
        elif testCase == "open":
            dispIncr = numpy.array([[ 9.1, 10.1,],
                                    [ 9.2, 10.2,],
                                    [ 9.3, 10.3,],
                                    [ 9.4, 10.4,],
                                    [ 9.5, 10.5,],
                                    [ 9.6, 10.6,],
                                    [ 9.7, 10.7,],
                                    [ 9.9, 10.9,],
                                    [ 9.8, 10.8,],
                                    [ 9.0, 10.0,],])

    # ------------------------------------------------------------------
    fieldTpdt = fieldT + fieldIncr

    tractionShear = abs(fieldTpdt[:,0]) / area
    tractionNormal = fieldTpdt[:,1] / area

    print "tractionShear",tractionShear
    print "tractionNormal",tractionNormal

    friction = -0.6 * tractionNormal;

    print "friction",friction

    lagrangeTpdt0 = friction * fieldTpdt[:,0] / tractionShear

    lagrangeIncr0 = lagrangeTpdt0 - fieldT[:,0]

    print "lagrangeIncr0",lagrangeIncr0

    dlagrange0 = (tractionShear - friction) * fieldTpdt[:,0] / tractionShear
    
    print "dlagrange0",dlagrange0

    if testCase == "slip": 
        dLagrange = numpy.vstack((dlagrange0, dlagrange1))
        dLagrange = numpy.transpose(dLagrange)
        dLagrange = numpy.reshape(dLagrange, m)
    elif testCase == "open":
        dLagrange = numpy.reshape(disp+dispIncr, n)
        dLagrange = dLagrange[indexL]

    print "dLagrange \n", dLagrange

    RHS = numpy.dot(numpy.transpose(C),dLagrange)
    duN = -numpy.dot(inv(jacobianN),RHS)
    duP = numpy.dot(inv(jacobianP),RHS)
    
    dispRel = duP - duN

    slipVertex = numpy.dot(C,dispRel)
    slipVertex = numpy.reshape(slipVertex, (m/DOF,DOF))
    if testCase == "slip":
        slipVertex[:,1] = 0
    mask = slipVertex[:,1] < 0.0
    slipVertex[mask,1] = 0
    slipVertex = numpy.reshape(slipVertex, m)

    print "duN \n", duN
    print "duP \n", duP

    dispIncrE = dispIncr
    dispIncrE = numpy.reshape(dispIncrE, n)

    dispIncrE[indexL] = dispIncrE[indexL] - dLagrange
    dispIncrE[indexN] = dispIncrE[indexN] - \
        0.5*numpy.dot(C.transpose(), slipVertex)
    dispIncrE[indexP] = dispIncrE[indexP] + \
        0.5*numpy.dot(C.transpose(), slipVertex)

    dispIncrE = numpy.reshape(dispIncrE, (n/DOF,DOF))
    slipVertex = numpy.reshape(slipVertex, (m/DOF,DOF))

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
        fieldIncr = numpy.array([[8.7, 9.7, -10.7],
                                 [8.9, 9.9, -10.9],
                                 [8.1, 9.1, -10.1]])
        area = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0])
        
        jacobianN = numpy.array(
            [[1.0, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7],
             [13.1, 1.0, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8], 
             [16.2, 16.3, 1.0, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9],
             [19.3, 19.4, 19.5, 1.0, 19.6, 19.7, 19.8, 19.9, 20.0],
             [22.4, 22.5, 22.6, 22.7, 1.0, 22.8, 22.9, 23.0, 23.1],
             [25.5, 25.6, 25.7, 25.8, 25.9, 1.0, 26.0, 26.1, 26.2],
             [28.6, 28.7, 28.8, 28.9, 29.0, 29.1, 1.0, 29.2, 29.3],
             [31.7, 31.8, 31.9, 32.0, 32.1, 32.2, 32.3, 1.0, 32.4],
             [34.8, 34.9, 35.0, 35.1, 35.2, 35.3, 35.4, 35.5, 1.0]])

        jacobianP = numpy.array(
            [[  1.0, 48.7, 48.8, 48.9, 49.0, 49.1, 49.2, 49.3, 49.4,],
             [ 51.8,  1.0, 51.9, 52.0, 52.1, 52.2, 52.3, 52.4, 52.5,],
             [ 54.9, 55.0,  1.0, 55.1, 55.2, 55.3, 55.4, 55.5, 55.6,],
             [ 58.0, 58.1, 58.2,  1.0, 58.3, 58.4, 58.5, 58.6, 58.7,],
             [ 61.1, 61.2, 61.3, 61.4,  1.0, 61.5, 61.6, 61.7, 61.8,],
             [ 64.2, 64.3, 64.4, 64.5, 64.6,  1.0, 64.7, 64.8, 64.9,],
             [ 67.3, 67.4, 67.5, 67.6, 67.7, 67.8,  1.0, 67.9, 68.0,],
             [ 70.4, 70.5, 70.6, 70.7, 70.8, 70.9, 71.0,  1.0, 71.1,],
             [ 73.5, 73.6, 73.7, 73.8, 73.9, 74.0, 74.1, 74.2,  1.0,],])

        disp = numpy.array([[ 7.1, 8.1, 9.1,],
                            [ 7.2, 8.2, 9.2,],
                            [ 7.3, 8.3, 9.3,],
                            [ 7.4, 8.4, 9.4,],
                            [ 7.5, 8.5, 9.5,],
                            [ 7.6, 8.6, 9.6,],
                            [ 7.8, 8.8, 9.8,],
                            [ 7.0, 8.0, 9.0,],
                            [ 7.7, 8.7, 9.7,],
                            [ 7.9, 8.9, 9.9,],
                            [ 7.1, 8.1, 9.1,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 8.1, 9.1, 10.1,],
                                    [ 8.2, 9.2, 10.2,],
                                    [ 8.3, 9.3, 10.3,],
                                    [ 8.4, 9.4, 10.4,],
                                    [ 8.5, 9.5, 10.5,],
                                    [ 8.6, 9.6, 10.6,],
                                    [ 8.8, 9.8, 10.8,],
                                    [ 8.0, 9.0, 10.0,],
                                    [ 8.7, 9.7, -10.7,],
                                    [ 8.9, 9.9, -10.9,],
                                    [ 8.1, 9.1, -10.1,],])            
        elif testCase == "open":
            dispIncr = numpy.array([[ 8.1, 9.1, 10.1,],
                                    [ 8.2, 9.2, 10.2,],
                                    [ 8.3, 9.3, 10.3,],
                                    [ 8.4, 9.4, 10.4,],
                                    [ 8.5, 9.5, 10.5,],
                                    [ 8.6, 9.6, 10.6,],
                                    [ 8.8, 9.8, 10.8,],
                                    [ 8.0, 9.0, 10.0,],
                                    [ 8.7, 9.7, 10.7,],
                                    [ 8.9, 9.9, 10.9,],
                                    [ 8.1, 9.1, 10.1,],])


    elif cell == "hex8":
        dlagrange2 = numpy.zeros(4)
        indexL = numpy.arange(48,60)
        indexN = numpy.arange(12,24)
        indexP = numpy.arange(36,48)
        n = 60
        m = 12
        DOF = 3

        fieldT = numpy.array([[5.4, 7.4, 9.4],
                              [5.6, 7.6, 9.6],
                              [5.8, 7.8, 9.8],
                              [5.0, 7.0, 9.0]])
        fieldIncr = numpy.array([[-6.4, -8.4, -10.4],
                                 [-6.6, -8.6, -10.6],
                                 [-6.8, -8.8, -10.8],
                                 [-6.0, -8.0, -10.0]])
        area = numpy.array([1.0, 1.0, 1.0, 1.0])
    
        jacobianN = numpy.array(
            [[   1.0,  72.1,  72.2,  72.3,  72.4,  72.5,  72.6,  72.7,  72.8,  72.9,  73.0,  73.1,],
             [  77.9,   1.0,  78.0,  78.1,  78.2,  78.3,  78.4,  78.5,  78.6,  78.7,  78.8,  78.9,],
             [  83.7,  83.8,   1.0,  83.9,  84.0,  84.1,  84.2,  84.3,  84.4,  84.5,  84.6,  84.7,],
             [  89.5,  89.6,  89.7,   1.0,  89.8,  89.9,  90.0,  90.1,  90.2,  90.3,  90.4,  90.5,],
             [  95.3,  95.4,  95.5,  95.6,   1.0,  95.7,  95.8,  95.9,  96.0,  96.1,  96.2,  96.3,],
             [ 101.1, 101.2, 101.3, 101.4, 101.5,   1.0, 101.6, 101.7, 101.8, 101.9, 102.0, 102.1,],
             [ 106.9, 107.0, 107.1, 107.2, 107.3, 107.4,   1.0, 107.5, 107.6, 107.7, 107.8, 107.9,],
             [ 112.7, 112.8, 112.9, 113.0, 113.1, 113.2, 113.3,   1.0, 113.4, 113.5, 113.6, 113.7,],
             [ 118.5, 118.6, 118.7, 118.8, 118.9, 119.0, 119.1, 119.2,   1.0, 119.3, 119.4, 119.5,],
             [ 124.3, 124.4, 124.5, 124.6, 124.7, 124.8, 124.9, 125.0, 125.1,   1.0, 125.2, 125.3,],
             [ 130.1, 130.2, 130.3, 130.4, 130.5, 130.6, 130.7, 130.8, 130.9, 131.0,   1.0, 131.1,],
             [ 135.9, 136.0, 136.1, 136.2, 136.3, 136.4, 136.5, 136.6, 136.7, 136.8, 136.9,   1.0,],])


        jacobianP = numpy.array(
            [[   1.0, 214.9, 215.0, 215.1, 215.2, 215.3, 215.4, 215.5, 215.6, 215.7, 215.8, 215.9,],
             [ 220.7,   1.0, 220.8, 220.9, 221.0, 221.1, 221.2, 221.3, 221.4, 221.5, 221.6, 221.7,],
             [ 226.5, 226.6,   1.0, 226.7, 226.8, 226.9, 227.0, 227.1, 227.2, 227.3, 227.4, 227.5,],
             [ 232.3, 232.4, 232.5,   1.0, 232.6, 232.7, 232.8, 232.9, 233.0, 233.1, 233.2, 233.3,],
             [ 238.1, 238.2, 238.3, 238.4,   1.0, 238.5, 238.6, 238.7, 238.8, 238.9, 239.0, 239.1,],
             [ 243.9, 244.0, 244.1, 244.2, 244.3,   1.0, 244.4, 244.5, 244.6, 244.7, 244.8, 244.9,],
             [ 249.7, 249.8, 249.9, 250.0, 250.1, 250.2,   1.0, 250.3, 250.4, 250.5, 250.6, 250.7,],
             [ 255.5, 255.6, 255.7, 255.8, 255.9, 256.0, 256.1,   1.0, 256.2, 256.3, 256.4, 256.5,],
             [ 261.3, 261.4, 261.5, 261.6, 261.7, 261.8, 261.9, 262.0,   1.0, 262.1, 262.2, 262.3,],
             [ 267.1, 267.2, 267.3, 267.4, 267.5, 267.6, 267.7, 267.8, 267.9,   1.0, 268.0, 268.1,],
             [ 272.9, 273.0, 273.1, 273.2, 273.3, 273.4, 273.5, 273.6, 273.7, 273.8,   1.0, 273.9,],
             [ 278.7, 278.8, 278.9, 279.0, 279.1, 279.2, 279.3, 279.4, 279.5, 279.6, 279.7,   1.0,],])

        disp = numpy.array([[ 4.1, 6.1, 8.1,],
                            [ 4.2, 6.2, 8.2,],
                            [ 4.3, 6.3, 8.3,],
                            [ 4.4, 6.4, 8.4,],
                            [ 4.5, 6.5, 8.5,],
                            [ 4.6, 6.6, 8.6,],
                            [ 4.7, 6.7, 8.7,],
                            [ 4.8, 6.8, 8.8,],
                            [ 4.9, 6.9, 8.9,],
                            [ 4.0, 6.0, 8.0,],
                            [ 5.1, 7.1, 9.1,],
                            [ 5.2, 7.2, 9.2,],
                            [ 5.3, 7.3, 9.3,],
                            [ 5.5, 7.5, 9.5,],
                            [ 5.7, 7.7, 9.7,],
                            [ 5.9, 7.9, 9.9,],
                            [ 5.4, 7.4, 9.4,],
                            [ 5.6, 7.6, 9.6,],
                            [ 5.8, 7.8, 9.8,],
                            [ 5.0, 7.0, 9.0,],])

        if testCase == "slip":
            dispIncr = numpy.array([[ 5.1, 7.1, 9.1,],
                                    [ 5.2, 7.2, 9.2,],
                                    [ 5.3, 7.3, 9.3,],
                                    [ 5.4, 7.4, 9.4,],
                                    [ 5.5, 7.5, 9.5,],
                                    [ 5.6, 7.6, 9.6,],
                                    [ 5.7, 7.7, 9.7,],
                                    [ 5.8, 7.8, 9.8,],
                                    [ 5.9, 7.9, 9.9,],
                                    [ 5.0, 7.0, 9.0,],
                                    [ 6.1, 8.1, 10.1,],
                                    [ 6.2, 8.2, 10.2,],
                                    [ 6.3, 8.3, 10.3,],
                                    [ 6.5, 8.5, 10.5,],
                                    [ 6.7, 8.7, 10.7,],
                                    [ 6.9, 8.9, 10.9,],
                                    [ -6.4, -8.4, -10.4,],
                                    [ -6.6, -8.6, -10.6,],
                                    [ -6.8, -8.8, -10.8,],
                                    [ -6.0, -8.0, -10.0,],])          
        elif testCase == "open":
            dispIncr = numpy.array([[ 5.1, 7.1, 9.1,],
                                    [ 5.2, 7.2, 9.2,],
                                    [ 5.3, 7.3, 9.3,],
                                    [ 5.4, 7.4, 9.4,],
                                    [ 5.5, 7.5, 9.5,],
                                    [ 5.6, 7.6, 9.6,],
                                    [ 5.7, 7.7, 9.7,],
                                    [ 5.8, 7.8, 9.8,],
                                    [ 5.9, 7.9, 9.9,],
                                    [ 5.0, 7.0, 9.0,],
                                    [ 6.1, 8.1, 10.1,],
                                    [ 6.2, 8.2, 10.2,],
                                    [ 6.3, 8.3, 10.3,],
                                    [ 6.5, 8.5, 10.5,],
                                    [ 6.7, 8.7, 10.7,],
                                    [ 6.9, 8.9, 10.9,],
                                    [ 6.4, 8.4, 10.4,],
                                    [ 6.6, 8.6, 10.6,],
                                    [ 6.8, 8.8, 10.8,],
                                    [ 6.0, 8.0, 10.0,],])

    # ------------------------------------------------------------------
    fieldTpdt = fieldT + fieldIncr

    tractionShear = (fieldTpdt[:,0]**2 + fieldTpdt[:,1]**2)**0.5 / area
    tractionNormal = fieldTpdt[:,2] / area

    print "tractionShear",tractionShear
    print "tractionNormal",tractionNormal

    friction = -0.6 * tractionNormal;

    print "friction",friction

    lagrangeTpdt0 = friction * fieldTpdt[:,0] / tractionShear
    lagrangeTpdt1 = friction * fieldTpdt[:,1] / tractionShear

    lagrangeIncr0 = lagrangeTpdt0 - fieldT[:,0]
    lagrangeIncr1 = lagrangeTpdt1 - fieldT[:,1]

    print "lagrangeIncr0",lagrangeIncr0
    print "lagrangeIncr1",lagrangeIncr1

    dlagrange0 = (tractionShear - friction) * fieldTpdt[:,0] / tractionShear
    dlagrange1 = (tractionShear - friction) * fieldTpdt[:,1] / tractionShear
    
    print "dlagrange0",dlagrange0
    print "dlagrange1",dlagrange1

    D = numpy.array([[ 0, -1, 0,],
                     [ 0, 0, +1,],
                     [ -1, 0, 0,],])

    Z = numpy.zeros([3,3])

    if cell == "tet4":
        C1 = numpy.hstack((D, Z, Z))
        C2 = numpy.hstack((Z, D, Z))
        C3 = numpy.hstack((Z, Z, D))
        C = numpy.vstack((C1, C2, C3))
    elif cell == "hex8":
        C1 = numpy.hstack((D, Z, Z, Z))
        C2 = numpy.hstack((Z, D, Z, Z))
        C3 = numpy.hstack((Z, Z, D, Z))
        C4 = numpy.hstack((Z, Z, Z, D))        
        C = numpy.vstack((C1, C2, C3, C4))

    if testCase == "slip": 
        dLagrange = numpy.vstack((dlagrange0, dlagrange1, dlagrange2))
        dLagrange = numpy.transpose(dLagrange)
        dLagrange = numpy.reshape(dLagrange, m)
    elif testCase == "open":
        dLagrange = numpy.reshape(disp+dispIncr, n)
        dLagrange = dLagrange[indexL]

    print "dLagrange \n", dLagrange

    RHS = numpy.dot(numpy.transpose(C),dLagrange)
    duN = -numpy.dot(inv(jacobianN),RHS)
    duP = numpy.dot(inv(jacobianP),RHS)
    
    dispRel = duP - duN

    slipVertex = numpy.dot(C,dispRel)
    slipVertex = numpy.reshape(slipVertex, (m/DOF,DOF))
    if testCase == "slip":
        slipVertex[:,2] = 0    
    mask = slipVertex[:,2] < 0.0
    slipVertex[mask,2] = 0
    slipVertex = numpy.reshape(slipVertex, m)

    print "duN \n", duN
    print "duP \n", duP

    dispIncrE = dispIncr
    dispIncrE = numpy.reshape(dispIncrE, n)

    dispIncrE[indexL] = dispIncrE[indexL] - dLagrange
    dispIncrE[indexN] = dispIncrE[indexN] - \
        0.5*numpy.dot(C.transpose(), slipVertex)
    dispIncrE[indexP] = dispIncrE[indexP] + \
        0.5*numpy.dot(C.transpose(), slipVertex)

    dispIncrE = numpy.reshape(dispIncrE, (n/DOF,DOF))
    slipVertex = numpy.reshape(slipVertex, (m/DOF,DOF))

    print "dispIncrE\n", printdata(dispIncrE)
    print "slipVertexE\n", printdata(slipVertex)

