cell = "hex8"
dim = "3d"

import numpy

# ----------------------------------------------------------------------
if dim == "2d":
    if cell == "tri3":

        fieldT = numpy.array([[8.6, 9.6],
                              [8.8, 9.8]])
        fieldIncr = numpy.array([[9.6, -10.6],
                                 [9.8, -10.8]])
        area = numpy.array([1.0, 1.0])
    

    elif cell == "quad4":
        fieldT = numpy.array([[8.8, 9.8],
                              [8.0, 9.0]])
        fieldIncr = numpy.array([[-9.8, -10.8],
                                 [-9.0, -10.0]])
        area = numpy.array([1.0, 1.0])
    

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

    slipVertex0 = 2 * dlagrange0

    print "slipVertex0",slipVertex0

    print "fieldIncr",numpy.array([lagrangeIncr0]).transpose()
    print "slip",numpy.array([slipVertex0]).transpose()

# ----------------------------------------------------------------------
elif dim == "3d":
    if cell == "tet4":

        fieldT = numpy.array([[7.7, 8.7, 9.7],
                              [7.9, 8.9, 9.9],
                              [7.1, 8.1, 9.1]])
        fieldIncr = numpy.array([[8.7, 9.7, -10.7],
                                 [8.9, 9.9, -10.9],
                                 [8.1, 9.1, -10.1]])
        area = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0])
    

    elif cell == "hex8":
        fieldT = numpy.array([[5.4, 7.4, 9.4],
                              [5.6, 7.6, 9.6],
                              [5.8, 7.8, 9.8],
                              [5.0, 7.0, 9.0]])
        fieldIncr = numpy.array([[-6.4, -8.4, -10.4],
                                 [-6.6, -8.6, -10.6],
                                 [-6.8, -8.8, -10.8],
                                 [-6.0, -8.0, -10.0]])
        area = numpy.array([1.0, 1.0, 1.0, 1.0])
    

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

    slipVertex0 = 2 * dlagrange0
    slipVertex1 = 2 * dlagrange1

    print "slipVertex0",slipVertex0
    print "slipVertex1",slipVertex1

    print "fieldIncr",numpy.array([lagrangeIncr0, lagrangeIncr1]).transpose()
    print "slip",numpy.array([slipVertex0, slipVertex1]).transpose()

