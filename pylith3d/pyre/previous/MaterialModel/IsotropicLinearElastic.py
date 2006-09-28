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
# $Id: IsotropicLinearElastic.py,v 1.1 2004/09/23 17:39:51 willic3 Exp $

# End of file 
