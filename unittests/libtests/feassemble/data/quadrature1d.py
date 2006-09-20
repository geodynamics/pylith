#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

def N0(p):
    return -0.5*p*(1-p)


def N1(p):
    return (1-p**2)


def N2(p):
    return 0.5*p*(1+p)


def N0p(p):
    return +0.5*p - 0.5*(1-p)


def N1p(p):
    return -2.0*p


def N2p(p):
    return 0.5*p + 0.5*(1+p)



q = [-1.0/3.0**0.5, 1.0/3.0**0.5]
#q = [-1.0, 0.0, 1.0]
print "basis = {"
for qp in q:
    print "%16.8e, %16.8e, %16.8e," % (N0(qp), N1(qp), N2(qp))
print "}"

print "basisDeriv = {"
for qp in q:
    print "%16.8e, %16.8e, %16.8e," % (N0p(qp), N1p(qp), N2p(qp))
print "}"


# End of file
