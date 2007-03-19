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

## @file pylith-0.8/pylith3d/pylith3d/MaterialModel/IsotropicLinearGenMaxwellViscoelastic.py

## @brief Python PyLith-0.8 generalized Maxwell Model.

from pylith3d.MaterialModel.MaterialModel import MaterialModel

# IsotropicLinearGenMaxwellViscoelastic class
class IsotropicLinearGenMaxwellViscoelastic(MaterialModel):
    """
    Python PyLith-0.8 material model consisting of 3 Maxwell elements in parallel.
    For now, the standard linear solid may be approximated by setting one viscosity to
    a large value and setting one of the shear ratios to zero.
    """
        
    def __init__(self):
        MaterialModel.__init__(self)

        import journal
        self.trace = journal.debug("pylith3d.trace")
        self.trace.log("Hello from IsotropicLinearGenMaxwellViscoelastic.__init__!")

        self.materialModel = 6
        self.numberProperties = 9
        # materialVariation flag is set to false since the material matrix does not vary
        # with the state variables.  It does vary with the time step size, however.
        self.materialVariation = False
        self.numberStateVariables = 30
        # This model assumes exactly 3 Maxwell models.
        self.propertyDict = {'density': None,
                             'youngsModulus': None,
                             'poissonsRatio': None,
                             'shearRatio1': None,
                             'viscosity1': None,
                             'shearRatio2': None,
                             'viscosity2': None,
                             'shearRatio3': None,
                             'viscosity3': None}
        self.propertyPosition = ['density',
                                 'youngsModulus',
                                 'poissonsRatio',
                                 'shearRatio1',
                                 'viscosity1',
                                 'shearRatio2',
                                 'viscosity2',
                                 'shearRatio3',
                                 'viscosity3']
        self.propertyList = [0.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0]
        return

# version
# $Id: IsotropicLinearGenMaxwellViscoelastic.py,v 1.2 2004/08/12 16:49:07 willic3 Exp $

# End of file 
