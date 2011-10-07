test = "hex8"
vertex = 3

# ----------------------------------------------------------------------
if test == "line2":
    jL = 1.0; rL = +1.69546413727;
    jN = 2.2; duN = 1.2
    jP = 2.4; duP = 1.4

    Sinv = 1.0/(jL**2 * (1.0/jN + 1.0/jP))
    duL = Sinv * (-rL + jL*(duP-duN));

    dduN = jL / jN * duL
    dduP = -jL / jP * duL

    #print "Aru",Aru
    #print "Aruslip",Aruslip
    #print "Si",Si
    #print "dlip",dlp

    print duN+dduN,duP+dduP,duL

# ----------------------------------------------------------------------
elif test == "tri3" or test == "quad4":

    if test == "tri3":

        if vertex == 0:
            # Lagrange vertex 8, vertex N: 3, vertex P: 6
            jL = 1.0
            rLx = -(8.5-8.2) - (0.08241148423)
            rLy = -(9.5-9.2) - (1.89546413727)
            jN = 1.2; duNx = 3.2; duNy = 4.2;
            jP = 1.5; duPx = 3.5; duPy = 4.5;

        elif vertex == 1:
            # Lagrange vertex 9, vertex N: 4, vertex P: 7
            jL = 1.0
            rLx = -(8.7-8.3) - (0.14794836271)
            rLy = -(9.7-9.3) - (1.77538035254)
            jN = 1.3; duNx = 3.3; duNy = 4.3;
            jP = 1.7; duPx = 3.7; duPy = 4.7;
            
    elif test == "quad4":

        if vertex == 0:
            # Lagrange vertex 10, vertex N: 4, vertex P: 8
            jL = 1.0
            rLx = -(8.7-8.3) + 0.14794836271
            rLy = -(9.7-9.3) + 1.77538035254
            jN = 1.3; duNx = 3.3; duNy = 4.3;
            jP = 1.7; duPx = 3.7; duPy = 4.7;
            
        elif vertex == 1:
            # Lagrange vertex 11, vertex N: 5, vertex P: 9
            jL = 1.0
            rLx = -(8.9-8.4) + 0.08241148423
            rLy = -(9.9-9.4) + 1.89546413727
            jN = 1.4; duNx = 3.4; duNy = 4.4;
            jP = 1.9; duPx = 3.9; duPy = 4.9;

    Sinv = 1.0/(jL**2 * (1.0/jN + 1.0/jP))
    duLx = Sinv * (-rLx + jL*(duPx - duNx))
    duLy = Sinv * (-rLy + jL*(duPy - duNy))

    dduNx = jL / jN * duLx
    dduNy = jL / jN * duLy

    dduPx = -jL / jP * duLx
    dduPy = -jL / jP * duLy
            
    print duNx+dduNx,duNy+dduNy
    print duPx+dduPx,duPy+dduPy
    print duLx,duLy

# ----------------------------------------------------------------------
elif test == "tet4" or test == "hex8":

    if test == "tet4":

        if vertex == 0:
            # Lagrange vertex 10, vertex N: 3, vertex P: 7
            jL = 1.0/3.0
            rLx = -1.0/3.0*(7.6-7.2 + -0.07938069066)
            rLy = -1.0/3.0*(8.6-8.2 + -1.82575588523)
            rLz = -1.0/3.0*(9.6-9.2 + 0.55566483464)
            jN = 1.2; duNx = 3.2; duNy = 4.2; duNz = 5.2;
            jP = 1.6; duPx = 3.6; duPy = 4.6; duPz = 5.6;

        elif vertex == 1:
            # Lagrange vertex 11, vertex N: 4, vertex P: 8
            jL = 1.0/3.0
            rLx = -1.0/3.0*(7.8-7.3 + -0.14140241667)
            rLy = -1.0/3.0*(8.8-8.3 + -1.69682900001)
            rLz = -1.0/3.0*(9.8-9.3 + 0.56560966667)
            jN = 1.3; duNx = 3.3; duNy = 4.3; duNz = 5.3;
            jP = 1.8; duPx = 3.8; duPy = 4.8; duPz = 5.8;
            
        elif vertex == 2:
            # Lagrange vertex 12, vertex N: 5, vertex P: 9
            jL = 1.0/3.0
            rLx = -1.0/3.0*(7.0-7.4 + -0.18205179147)
            rLy = -1.0/3.0*(8.0-8.4 + -1.51709826228)
            rLz = -1.0/3.0*(9.0-9.4 + 0.54615537442)
            jN = 1.4; duNx = 3.4; duNy = 4.4; duNz = 5.4;
            jP = 1.0; duPx = 3.0; duPy = 4.0; duPz = 5.0;
            
    elif test == "hex8":

        if vertex == 0:
            # Lagrange vertex 18, vertex N: 6, vertex P: 14
            jL = 1.0
            rLx = -(5.3-4.5+0.07938069066)
            rLy = -(7.3-6.5+1.82575588523)
            rLz = -(9.3-8.5+0.55566483464)
            jN = 1.5; duNx = 3.5; duNy = 4.5; duNz = 5.5;
            jP = 1.3; duPx = 3.3; duPy = 4.3; duPz = 5.3;
            
        elif vertex == 1:
            # Lagrange vertex 19, vertex N: 7, vertex P: 15
            jL = 1.0
            rLx = -(5.5-4.6+0.14140241667)
            rLy = -(7.5-6.6+1.69682900001)
            rLz = -(9.5-8.6+0.56560966667)
            jN = 1.6; duNx = 3.6; duNy = 4.6; duNz = 5.6;
            jP = 1.5; duPx = 3.5; duPy = 4.5; duPz = 5.5;
            
        elif vertex == 2:
            # Lagrange vertex 20, vertex N: 8, vertex P: 16
            jL = 1.0
            rLx = -(5.7-4.7+0.18205179147)
            rLy = -(7.7-6.7+1.51709826228)
            rLz = -(9.7-8.7+0.54615537442)
            jN = 1.7; duNx = 3.7; duNy = 4.7; duNz = 5.7;
            jP = 1.7; duPx = 3.7; duPy = 4.7; duPz = 5.7;
            
        elif vertex == 3:
            # Lagrange vertex 21, vertex N: 9, vertex P: 17
            jL = 1.0
            rLx = -(5.9-4.8+0.19904410828)
            rLy = -(7.9-6.8+1.29378670385)
            rLz = -(9.9-8.8+0.49761027071)
            jN = 1.8; duNx = 3.8; duNy = 4.8; duNz = 5.8;
            jP = 1.9; duPx = 3.9; duPy = 4.9; duPz = 5.9;
            
    Sinv = 1.0/(jL**2 * (1.0/jN + 1.0/jP))
    duLx = Sinv * (-rLx + jL*(duPx - duNx))
    duLy = Sinv * (-rLy + jL*(duPy - duNy))
    duLz = Sinv * (-rLz + jL*(duPz - duNz))

    dduNx = jL / jN * duLx
    dduNy = jL / jN * duLy
    dduNz = jL / jN * duLz

    dduPx = -jL / jP * duLx
    dduPy = -jL / jP * duLy
    dduPz = -jL / jP * duLz
            
    print duNx+dduNx,duNy+dduNy, duNz+dduNz
    print duPx+dduPx,duPy+dduPy, duPz+dduPz
    print duLx,duLy, duLz

