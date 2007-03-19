#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

## @file pylith-0.8/pylith3d/pylith3d/MaterialModel/IsotropicPowerLawMaxwellViscoelastic.py

## @brief Python PyLith-0.8 power-law Maxwell viscoelastic model.

from pylith3d.MaterialModel.MaterialModel import MaterialModel

class IsotropicPowerLawMaxwellViscoelastic(MaterialModel):
    """
    Python PyLith-0.8 definitions for a power-law Maxwell viscoelastic material.
    This version uses the Bathe Effective Stress Function (ESF) approach.
    """

    def __init__(self):
        MaterialModel.__init__(self)

        import journal
        self.trace = journal.debug("pylith3d.trace")
        self.trace.log("Hello from IsotropicPowerLawMaxwellViscoelastic.__init__!")

        self.materialModel = 7
        self.numberProperties = 5
        self.materialVariation = True
        self.numberStateVariables = 18
        self.propertyDict = {'density': None,
                             'youngsModulus': None,
                             'poissonsRatio': None,
                             'powerLawExponent': None,
                             'viscosityCoefficient': None}
        self.propertyPosition = ['density',
                                 'youngsModulus',
                                 'poissonsRatio',
                                 'powerLawExponent',
                                 'viscosityCoefficient']
        self.propertyList = [0.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0]
        return

# version
# $Id: IsotropicPowerLawMaxwellViscoelastic.py,v 1.2 2004/08/12 16:49:07 willic3 Exp $

# End of file 
