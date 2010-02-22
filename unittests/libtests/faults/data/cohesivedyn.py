cell = "hex8"
dim = "3d"

import numpy

# ----------------------------------------------------------------------
if dim == "3d":
    if cell == "tet4":

        fieldT = numpy.array([[7.7, 8.7, 9.7],
                              [7.9, 8.9, 9.9],
                              [7.1, 8.1, 9.1]])
        fieldIncr = numpy.array([[8.7, 9.7, -10.7],
                                 [8.9, 9.9, -10.9],
                                 [8.1, 9.1, -10.1]])
        initialTract = numpy.array([[1.1, 2.1, -3.1],
                                    [1.1, 2.1, -3.1],
                                    [1.1, 2.1, -3.1]])
        area = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0])
    

    elif cell == "hex8":
        fieldT = numpy.array([[5.4, 7.4, 9.4],
                              [5.6, 7.6, 9.6],
                              [5.8, 7.8, 9.8],
                              [5.0, 7.0, 9.0]])
        fieldIncr = numpy.array([[6.4, 8.4, -10.4],
                                 [6.6, 8.6, -10.6],
                                 [6.8, 8.8, -10.8],
                                 [6.0, 8.0, -10.0]])
        initialTract = numpy.array([[1.063397471, 2.063397471, -3.063397471], 
                                    [1.121132498, 2.121132498, -3.121132498], 
                                    [1.178867525, 2.178867525, -3.178867525],
                                    [1.236602552, 2.236602552, -3.236602552]])
        area = numpy.array([1.0, 1.0, 1.0, 1.0])
    

        tp = numpy.array([9.4, 9.6, 9.8, 9])
        tdt = numpy.array([-10.4, -10.6, -10.8, -10])
        it = numpy.array([-3.0634, -3.12113, -3.17887, -3.2366])
        
        tps1 = numpy.array([5.4, 5.6, 5.8, 5])
        tps2 = numpy.array([7.4, 7.6, 7.8, 7])
        tdts1 = numpy.array([6.4, 6.6, 6.8, 6])
        tdts2 = numpy.array([8.4, 8.6, 8.8, 8])
        its1 = numpy.array([1.0634, 1.12113, 1.17887, 1.2366])
        its2 = numpy.array([2.0634, 2.12113, 2.17887, 2.2366])
        
        tpdts1 = tps1 + tdts1
        tpdts2 = tps2 + tdts2
        
        ts = (tpdts1**2 + tpdts2**2)**0.5
        nt = tp + tdt + it
        
        friction = -0.6 * nt
        
        dL0 = (ts - friction) * tpdts1 / ts;
        dL1 = (ts - friction) * tpdts2 / ts;
        
        slipVertex0 = 2 * dL0;
        slipVertex1 = 2 * dL1;

        tpdts1 = friction * tpdts1 / ts
        tpdts2 = friction * tpdts2 / ts
    

    # ------------------------------------------------------------------
    fieldTpdt = fieldT + fieldIncr

    tractionShear = (fieldTpdt[:,0]**2 + fieldTpdt[:,1]**2)**0.5 / area
    tractionNormal = initialTract[:,2] + fieldTpdt[:,2] / area

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

