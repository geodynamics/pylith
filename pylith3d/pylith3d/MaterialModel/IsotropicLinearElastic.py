#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2004  All Rights Reserved
#
#  Copyright 2004 Rensselaer Polytechnic Institute.
#  All worldwide rights reserved.  A license to use, copy, modify and
#  distribute this software for non-commercial research purposes only
#  is hereby granted, provided that this copyright notice and
#  accompanying disclaimer is not modified or removed from the software.
#
#  DISCLAIMER:  The software is distributed "AS IS" without any express
#  or implied warranty, including but not limited to, any implied
#  warranties of merchantability or fitness for a particular purpose
#  or any warranty of non-infringement of any current or pending patent
#  rights.  The authors of the software make no representations about
#  the suitability of this software for any particular purpose.  The
#  entire risk as to the quality and performance of the software is with
#  the user.  Should the software prove defective, the user assumes the
#  cost of all necessary servicing, repair or correction.  In
#  particular, neither Rensselaer Polytechnic Institute, nor the authors
#  of the software are liable for any indirect, special, consequential,
#  or incidental damages related to the software, to the maximum extent
#  the law permits.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pylith3d.MaterialModel.MaterialModel import MaterialModel

class IsotropicLinearElastic(MaterialModel):
    """Basic definitions for an isotropic linear elastic material"""

    def __init__(self):
        """Initialization for isotropic linear elastic material

        This corresponds to material model 1.

        There are 3 material properties:
            density
            Young's modulus
            Poisson's ratio

        The variable materialVariation is set to False, indicating that
        the material matrix does not vary with the state variables.  This
        means that only a single material matrix is required for each
        material group of this type.

        Only 2 state variables (stress and strain) are required to describe
        this material."""

        # print "Hello from IsotropicLinearElastic.__init__!"
        # print ""
        self.materialModel = 1
        self.numberProperties = 3
        self.materialVariation = False
        self.numberStateVariables = 2
        self.propertyDict = {'density': None,
                             'youngsModulus': None,
                             'poissonsRatio': None}
        self.propertyPosition = ['density',
                                 'youngsModulus',
                                 'poissonsRatio']
        self.propertyList = [0.0,
                             0.0,
                             0.0]
        return

# version
# $Id: IsotropicLinearElastic.py,v 1.2 2004/08/12 16:47:30 willic3 Exp $

# End of file 
