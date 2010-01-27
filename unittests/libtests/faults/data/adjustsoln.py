# Line2
C = 1.0
dk = 1.89546413727; dl = 1.5
Ai = 2.2; ri = 7.5; ui = 7.2; dui = 1.2
Aj = 2.4; rj = -7.5; uj = 7.4; duj = 1.4

Si = (Ai * Aj) / (Ai + Aj)
ddl = Si * (C*(ri/Ai - rj/Aj + ui - uj) - dk)

ddui = -C / Ai * ddl
dduj = +C / Aj * ddl

print ddl
print ddui
print dui+ddui,duj+dduj,dl+ddl
