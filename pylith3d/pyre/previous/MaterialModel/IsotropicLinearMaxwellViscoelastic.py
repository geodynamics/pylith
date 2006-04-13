#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pylith3d.MaterialModel.MaterialModel import MaterialModel

class IsotropicLinearMaxwellViscoelastic(MaterialModel):

    def __init__(self):
        # print "Hello from IsotropicLinearMaxwellViscoelastic.__init__!"
        # print ""
        self.materialModel = 5
        self.numberProperties = 4
        self.propertyDict = {'density': None,
                             'youngsModulus': None,
                             'poissonsRatio': None,
                             'viscosity': None}
        self.propertyPosition = ['density',
                                 'youngsModulus',
                                 'poissonsRatio',
                                 'viscosity']
        self.propertyList = [0.0,
                             0.0,
                             0.0,
                             0.0]
        return

# version
# $Id: IsotropicLinearMaxwellViscoelastic.py,v 1.1 2004/09/23 17:39:51 willic3 Exp $

# End of file 
