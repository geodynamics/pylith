test = "hex8"
vertex = 3

# ----------------------------------------------------------------------
if test == "line2":
    C = 1.0
    dk = 1.89546413727; l = 1.5
    Ai = 2.2; ri = 7.5; ui = 7.2; dui = 1.2
    Aj = 2.4; rj = -7.5; uj = 7.4; duj = 1.4

    Si = (Ai * Aj) / (Ai + Aj)
    Aru = ri/Ai - rj/Aj + ui - uj
    Aruslip = -C*Aru - dk 
    dlp = Si * Aruslip

    ddui = +C / Ai * dlp
    dduj = -C / Aj * dlp

    #print "Aru",Aru
    #print "Aruslip",Aruslip
    #print "Si",Si
    #print "dlip",dlp

    print dui+ddui,duj+dduj,l+dlp

# ----------------------------------------------------------------------
elif test == "tri3" or test == "quad4":

    if test == "tri3":

        Cpx = 0.0
        Cpy = -1.0
        Cqx = -1.0
        Cqy = 0.0

        if vertex == 0:
            # Lagrange vertex 8, vertex i: 3, vertex j: 6
            dkp = 1.89546413727; lp = 3.6
            dkq = 0.08241148423; lq = 4.6
            
            # vertex i
            Aix = 1.2; Aiy = 1.2
            rix = -9.6; riy = -8.6
            uix = 8.2; uiy = 9.2
            duix = 3.2; duiy = 4.2
            
            # vertex j
            Ajx = 1.5; Ajy = 1.5
            rjx = +9.6; rjy = +8.6
            ujx = 8.5; ujy = 9.5
            dujx = 3.5; dujy = 4.5
            
        elif vertex == 1:
            # Lagrange vertex 9, vertex i: 4, vertex j: 7
            dkp = 1.77538035254; lp = 3.8
            dkq = 0.14794836271; lq = 4.8
            
            # vertex i
            Aix = 1.3; Aiy = 1.3
            rix = -9.8; riy = -8.8
            uix = 8.3; uiy = 9.3
            duix = 3.3; duiy = 4.3
            
            # vertex j
            Ajx = 1.7; Ajy = 1.7
            rjx = +9.8; rjy = +8.8
            ujx = 8.7; ujy = 9.7
            dujx = 3.7; dujy = 4.7
            
    elif test == "quad4":

        Cpx = 0.0
        Cpy = 1.0
        Cqx = 1.0
        Cqy = 0.0
        
        if vertex == 0:
            # Lagrange vertex 10, vertex i: 4, vertex j: 8
            dkp = 1.77538035254; lp = 3.8
            dkq = 0.14794836271; lq = 4.8
            
            # vertex i: 4
            Aix = 1.3; Aiy = 1.3
            rix = +9.8; riy = +8.8
            uix = 8.3; uiy = 9.3
            duix = 3.3; duiy = 4.3
            
            # vertex j: 8
            Ajx = 1.7; Ajy = 1.7
            rjx = -9.8; rjy = -8.8
            ujx = 8.7; ujy = 9.7
            dujx = 3.7; dujy = 4.7
            
        elif vertex == 1:
            # Lagrange vertex 11, vertex i: 5, vertex j: 9
            dkp = 1.89546413727; lp = 3.0
            dkq = 0.08241148423; lq = 4.0
            
            # vertex i: 5
            Aix = 1.4; Aiy = 1.4
            rix = +9.0; riy = +8.0
            uix = 8.4; uiy = 9.4
            duix = 3.4; duiy = 4.4
            
            # vertex j: 9
            Ajx = 1.9; Ajy = 1.9
            rjx = -9.0; rjy = -8.0
            ujx = 8.9; ujy = 9.9
            dujx = 3.9; dujy = 4.9

    Sppi = Aix*Ajx / (Aix + Ajx)
    Sqqi = Aix*Ajx / (Aix + Ajx)

    Arux = rix / Aix - rjx / Ajx + uix - ujx
    Aruy = riy / Aiy - rjy / Ajy + uiy - ujy

    Arup = Cpx*Arux + Cpy*Aruy
    Aruq = Cqx*Arux + Cqy*Aruy
    Arupslip = -Arup - dkp 
    Aruqslip = -Aruq - dkq 

    dlp = Sppi * Arupslip
    dlq = Sqqi * Aruqslip

    dduix = +1.0/Aix * (Cpx*dlp + Cqx*dlq)
    dduiy = +1.0/Aiy * (Cpy*dlp + Cqy*dlq)

    ddujx = -1.0/Ajx * (Cpx*dlp + Cqx*dlq)
    ddujy = -1.0/Ajy * (Cpy*dlp + Cqy*dlq)

    print "Sppi",Sppi
    print "Sqqi",Sqqi

    print "Arup",Arup
    print "Aruq",Aruq

    print "Arupslip",Arupslip
    print "Aruqslip",Aruqslip

    print "dlp",dlp
    print "dlq",dlq

    print "dduix:",dduix
    print "dduiy:",dduiy

    print "ddujx:",ddujx
    print "ddujy:",ddujy

    print duix+dduix,duiy+dduiy
    print dujx+ddujx,dujy+ddujy
    print lp+dlp,lq+dlq

# ----------------------------------------------------------------------
elif test == "tet4" or test == "hex8":

    if test == "tet4":

        Cpx = 0.0
        Cpy = 1.0
        Cpz = 0.0
        Cqx = 0.0
        Cqy = 0.0
        Cqz = 1.0
        Crx = 1.0
        Cry = 0.0
        Crz = 0.0

        if vertex == 0:
            # Lagrange vertex 10, vertex i: 3, vertex j: 7
            dkp = 1.82575588523; lp = 3.7
            dkq = -0.55566483464; lq = 4.7
            dkr = 0.07938069066; lr = 5.7
            
            # vertex i: 3
            Aix = 1.2; Aiy = 1.2; Aiz = 1.2
            rix = +9.7; riy = +7.7; riz = +8.7
            uix = 7.2; uiy = 8.2; uiz = 9.2
            duix = 3.2; duiy = 4.2; duiz = 5.2
            
            # vertex j: 7
            Ajx = 1.6; Ajy = 1.6; Ajz = 1.6
            rjx = -9.7; rjy = -7.7; rjz = -8.7
            ujx = 7.6; ujy = 8.6; ujz = 9.6
            dujx = 3.6; dujy = 4.6; dujz = 5.6
            
        elif vertex == 1:
            # Lagrange vertex 11, vertex i: 4, vertex j: 8
            dkp = 1.69682900001; lp = 3.9
            dkq = -0.56560966667; lq = 4.9
            dkr = 0.14140241667; lr = 5.9
            
            # vertex i: 4
            Aix = 1.3; Aiy = 1.3; Aiz = 1.3
            rix = +9.9; riy = +7.9; riz = +8.9
            uix = 7.3; uiy = 8.3; uiz = 9.3
            duix = 3.3; duiy = 4.3; duiz = 5.3
            
            # vertex j: 8
            Ajx = 1.8; Ajy = 1.8; Ajz = 1.8
            rjx = -9.9; rjy = -7.9; rjz = -8.9
            ujx = 7.8; ujy = 8.8; ujz = 9.8
            dujx = 3.8; dujy = 4.8; dujz = 5.8
            
        elif vertex == 2:
            # Lagrange vertex 12, vertex i: 5, vertex j: 9
            dkp = 1.51709826228; lp = 3.1
            dkq = -0.54615537442; lq = 4.1
            dkr = 0.18205179147; lr = 5.1
            
            # vertex i: 5
            Aix = 1.4; Aiy = 1.4; Aiz = 1.4
            rix = +9.1; riy = +7.1; riz = +8.1
            uix = 7.4; uiy = 8.4; uiz = 9.4
            duix = 3.4; duiy = 4.4; duiz = 5.4
            
            # vertex j: 9
            Ajx = 1.0; Ajy = 1.0; Ajz = 1.0
            rjx = -9.1; rjy = -7.1; rjz = -8.1
            ujx = 7.0; ujy = 8.0; ujz = 9.0
            dujx = 3.0; dujy = 4.0; dujz = 5.0
            
    elif test == "hex8":

        Cpx = 0.0
        Cpy = -1.0
        Cpz = 0.0
        Cqx = 0.0
        Cqy = 0.0
        Cqz = -1.0
        Crx = -1.0
        Cry = 0.0
        Crz = 0.0

        if vertex == 0:
            # Lagrange vertex 18, vertex i: 6, vertex j: 14
            dkp = 1.82575588523; lp = 3.4
            dkq = -0.55566483464; lq = 4.4
            dkr = 0.07938069066; lr = 5.4
            
            # vertex i: 6
            Aix = 1.5; Aiy = 1.5; Aiz = 1.5
            rix = -9.4; riy = -5.4; riz = -7.4
            uix = 4.5; uiy = 6.5; uiz = 8.5
            duix = 3.5; duiy = 4.5; duiz = 5.5
            
            # vertex j: 14
            Ajx = 1.3; Ajy = 1.3; Ajz = 1.3
            rjx = +9.4; rjy = +5.4; rjz = +7.4
            ujx = 5.3; ujy = 7.3; ujz = 9.3
            dujx = 3.3; dujy = 4.3; dujz = 5.3
            
        elif vertex == 1:
            # Lagrange vertex 19, vertex i: 7, vertex j: 15
            dkp = 1.69682900001; lp = 3.6
            dkq = -0.56560966667; lq = 4.6
            dkr = 0.14140241667; lr = 5.6
            
            # vertex i: 7
            Aix = 1.6; Aiy = 1.6; Aiz = 1.6
            rix = -9.6; riy = -5.6; riz = -7.6
            uix = 4.6; uiy = 6.6; uiz = 8.6
            duix = 3.6; duiy = 4.6; duiz = 5.6
            
            # vertex j: 15
            Ajx = 1.5; Ajy = 1.5; Ajz = 1.5
            rjx = +9.6; rjy = +5.6; rjz = +7.6
            ujx = 5.5; ujy = 7.5; ujz = 9.5
            dujx = 3.5; dujy = 4.5; dujz = 5.5
            
        elif vertex == 2:
            # Lagrange vertex 20, vertex i: 8, vertex j: 16
            dkp = 1.51709826228; lp = 3.8
            dkq = -0.54615537442; lq = 4.8
            dkr = 0.18205179147; lr = 5.8
            
            # vertex i: 8
            Aix = 1.7; Aiy = 1.7; Aiz = 1.7
            rix = -9.8; riy = -5.8; riz = -7.8
            uix = 4.7; uiy = 6.7; uiz = 8.7
            duix = 3.7; duiy = 4.7; duiz = 5.7
            
            # vertex j: 16
            Ajx = 1.7; Ajy = 1.7; Ajz = 1.7
            rjx = +9.8; rjy = +5.8; rjz = +7.8
            ujx = 5.7; ujy = 7.7; ujz = 9.7
            dujx = 3.7; dujy = 4.7; dujz = 5.7
            
        elif vertex == 3:
            # Lagrange vertex 21, vertex i: 9, vertex j: 17
            dkp = 1.29378670385; lp = 3.0
            dkq = -0.49761027071; lq = 4.0
            dkr = 0.19904410828; lr = 5.0
            
            # vertex i: 9
            Aix = 1.8; Aiy = 1.8; Aiz = 1.8
            rix = -9.0; riy = -5.0; riz = -7.0
            uix = 4.8; uiy = 6.8; uiz = 8.8
            duix = 3.8; duiy = 4.8; duiz = 5.8
            
            # vertex j: 17
            Ajx = 1.9; Ajy = 1.9; Ajz = 1.9
            rjx = +9.0; rjy = +5.0; rjz = +7.0
            ujx = 5.9; ujy = 7.9; ujz = 9.9
            dujx = 3.9; dujy = 4.9; dujz = 5.9
            
    Sppi = Aix*Ajx / (Aix + Ajx)
    Sqqi = Aiy*Ajy / (Aiy + Ajy)
    Srri = Aiz*Ajz / (Aiz + Ajz)

    Arux = rix / Aix - rjx / Ajx + uix - ujx
    Aruy = riy / Aiy - rjy / Ajy + uiy - ujy
    Aruz = riz / Aiz - rjz / Ajz + uiz - ujz

    Arup = Cpx*Arux + Cpy*Aruy + Cpz*Aruz
    Aruq = Cqx*Arux + Cqy*Aruy + Cqz*Aruz
    Arur = Crx*Arux + Cry*Aruy + Crz*Aruz
    Arupslip = -Arup - dkp
    Aruqslip = -Aruq - dkq
    Arurslip = -Arur - dkr

    dlp = Sppi * Arupslip
    dlq = Sqqi * Aruqslip
    dlr = Srri * Arurslip

    dduix = +1.0/Aix * (Cpx*dlp + Cqx*dlq + Crx*dlr)
    dduiy = +1.0/Aiy * (Cpy*dlp + Cqy*dlq + Cry*dlr)
    dduiz = +1.0/Aiy * (Cpz*dlp + Cqz*dlq + Crz*dlr)

    ddujx = -1.0/Ajx * (Cpx*dlp + Cqx*dlq + Crx*dlr)
    ddujy = -1.0/Ajy * (Cpy*dlp + Cqy*dlq + Cry*dlr)
    ddujz = -1.0/Ajz * (Cpz*dlp + Cqz*dlq + Crz*dlr)

    print "Sppi",Sppi
    print "Sqqi",Sqqi
    print "Srri",Srri

    print "Arup",Arup
    print "Aruq",Aruq
    print "Arur",Arur

    print "Arupslip",Arupslip
    print "Aruqslip",Aruqslip
    print "Arurslip",Arurslip

    print "dlp",dlp
    print "dlq",dlq
    print "dlr",dlr

    print "dduix:",dduix
    print "dduiy:",dduiy
    print "dduiz:",dduiz

    print "ddujx:",ddujx
    print "ddujy:",ddujy
    print "ddujz:",ddujz

    print duix+dduix,duiy+dduiy,duiz+dduiz
    print dujx+ddujx,dujy+ddujy,dujz+ddujz
    print lp+dlp,lq+dlq,lr+dlr

