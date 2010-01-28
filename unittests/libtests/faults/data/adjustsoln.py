test = "tri3"

# ----------------------------------------------------------------------
if test == "line2":
    C = 1.0
    dk = 1.89546413727; l = 1.5
    Ai = 2.2; ri = 7.5; ui = 7.2; dui = 1.2
    Aj = 2.4; rj = -7.5; uj = 7.4; duj = 1.4

    Si = (Ai * Aj) / (Ai + Aj)
    Aru = ri/Ai - rj/Aj + ui - uj
    Aruslip = C*Aru - dk
    dlp = Si * Aruslip

    ddui = -C / Ai * dlp
    dduj = +C / Aj * dlp

    #print "Aru",Aru
    #print "Aruslip",Aruslip
    #print "Si",Si
    #print "dlip",dlp

    print dui+ddui,duj+dduj,l+dlp

# ----------------------------------------------------------------------
elif test == "tri3":

    # Lagrange vertex 8, vertex i: 3, vertex j: 6
    Cpx = 0.0
    Cpy = -1.0
    Cqx = -1.0
    Cqy = 0.0

    dkp = 1.89546413727; lp = 1.5
    dkq = 0.08241148423; lq = 1.5
    
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

    Sppi = Aix*Ajx / (Aix + Ajx)
    Sqqi = Aix*Ajx / (Aix + Ajx)

    Arux = rix / Aix - rjx / Ajx + uix - ujx
    Aruy = riy / Aiy - rjy / Ajy + uiy - ujy

    Arup = Cpx*Arux + Cpy*Aruy
    Aruq = Cqx*Arux + Cqy*Aruy
    Arupslip = Arup - dkp
    Aruqslip = Aruq - dkq

    dlp = Sppi * Arupslip
    dlq = Sqqi * Aruqslip

    dduix = -1.0/Aix * (Cpx*dlp + Cqx*dlq)
    dduiy = -1.0/Aiy * (Cpy*dlp + Cqy*dlq)

    ddujx = +1.0/Ajx * (Cpx*dlp + Cqx*dlq)
    ddujy = +1.0/Ajy * (Cpy*dlp + Cqy*dlq)

    #print "Aru",Aru
    #print "Aruslip",Aruslip
    #print "Si",Si
    #print "dlip",dlp

    print duix+dduix,duiy+dduiy
    print dujx+ddujx,dujy+ddujy
    print lp+dlp,lq+dlq
