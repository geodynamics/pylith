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
    Python PyLith-0.8 material model consisting of a linear elastic spring in parallel
    with N Maxwell elements.
    """

    def readprop(self,file):
        """
        We need to override the generic property reader since this
        one has to increase the sizes of dictionaries and lists
        depending on the input values.
        """

        from pylith3d.KeywordValueParse import KeywordValueParse
        lineparse = KeywordValueParse()

        # print "Hellof from IsotropicLinearGenMaxwellViscoelastic.readprop!"

        endMaterial = False
        count = 0
        while not endMaterial:
            line=file.readline()
            if not line: break
            keyvals = lineparse.parseline(line)
            if count == 0:
                if keyvals[0] != "numberMaxwellElements":
                    raise ValueError, "First entry for this material type must be numberMaxwellElements"
                else:
                    numberModels = keyvals[1]
                    self.numberProperties += 2*numberModels
                    newModels = 1
                    while newModels < numberModels + 1:
                        shear = "shearModulus" + str(newModels)
                        viscosity = "viscosity" + str(newModels)
                        self.propertyDict += {shear: None, viscosity: None}
                        self.propertyPosition += [shear, viscosity]
                        self.propertyList += [0.0, 0.0]
                        newModels += 1
            count += 1
            if keyvals[3]:
                exec keyvals[0] + '=' + `keyvals[1]`
        for propertyIndex in range(len(self.propertyPosition)):
            self.propertyDict[self.propertyPosition[propertyIndex]] = eval(self.propertyPosition[propertyIndex])
            self.propertyList[propertyIndex] = self.propertyDict[self.propertyPosition[propertyIndex]]

        test = None in self.propertyDict.values()
        if test:
            position = self.propertyDict.values().index(None)
            noneKey = self.propertyDict.keys()[position]
            raise ValueError, "No value assigned for property: %r" % noneKey

        return
                
        
    def __init__(self):
        # print "Hello from IsotropicLinearGenMaxwellViscoelastic.__init__!"
        # print ""
        self.materialModel = 7
        self.numberProperties = 5
        # This model starts with no Maxwell elements defined.
        # They are added as the material properties file is read.
        self.propertyDict = {'density': None,
                             'numberMaxwellElements': None,
                             'bulkModulus': None,
                             'shearModulus': None,
                             'shearModulusInfinity': None}
        self.propertyPosition = ['density',
                                 'numberMaxwellElements',
                                 'bulkModulus',
                                 'shearModulus',
                                 'shearModulusInfinity']
        self.propertyList = [0.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0]
        return

# version
# $Id: IsotropicLinearGenMaxwellViscoelastic.py,v 1.2 2004/08/12 16:49:07 willic3 Exp $

# End of file 
