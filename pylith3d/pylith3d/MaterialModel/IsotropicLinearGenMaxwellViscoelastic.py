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
    # """
    # Python PyLith-0.8 material model consisting of a linear elastic spring in parallel
    # with N Maxwell elements.
    """

    # To prevent having to rewrite a bunch of other code sections, I am temporarily commenting
    # out the specialized readprop function and setting the number of Maxwell elements to 3.
    #temp def readprop(self,file):
    #temp     """
    #temp     We need to override the generic property reader since this
    #temp     one has to increase the sizes of dictionaries and lists
    #temp     depending on the input values.
    #temp     """

    #temp     from pylith3d.KeywordValueParse import KeywordValueParse
    #temp     lineparse = KeywordValueParse()

    #temp     # print "Hellof from IsotropicLinearGenMaxwellViscoelastic.readprop!"

    #temp     endMaterial = False
    #temp     count = 0
    #temp     while not endMaterial:
    #temp         line=file.readline()
    #temp         if not line: break
    #temp         keyvals = lineparse.parseline(line)
    #temp         if count == 0:
    #temp             if keyvals[0] != "numberMaxwellElements":
    #temp                 raise ValueError, "First entry for this material type must be numberMaxwellElements"
    #temp             else:
    #temp                 numberModels = keyvals[1]
    #temp                 self.numberProperties += 2*numberModels
    #temp                 newModels = 1
    #temp                 while newModels < numberModels + 1:
    #temp                     shear = "shearModulus" + str(newModels)
    #temp                     viscosity = "viscosity" + str(newModels)
    #temp                     self.propertyDict += {shear: None, viscosity: None}
    #temp                     self.propertyPosition += [shear, viscosity]
    #temp                     self.propertyList += [0.0, 0.0]
    #temp                     newModels += 1
    #temp         count += 1
    #temp         if keyvals[3]:
    #temp             exec keyvals[0] + '=' + `keyvals[1]`
    #temp     for propertyIndex in range(len(self.propertyPosition)):
    #temp         self.propertyDict[self.propertyPosition[propertyIndex]] = eval(self.propertyPosition[propertyIndex])
    #temp         self.propertyList[propertyIndex] = self.propertyDict[self.propertyPosition[propertyIndex]]

    #temp     test = None in self.propertyDict.values()
    #temp     if test:
    #temp         position = self.propertyDict.values().index(None)
    #temp         noneKey = self.propertyDict.keys()[position]
    #temp         raise ValueError, "No value assigned for property: %r" % noneKey

    #temp     return
                
        
    def __init__(self):
        # print "Hello from IsotropicLinearGenMaxwellViscoelastic.__init__!"
        # print ""
        self.materialModel = 7
        self.numberProperties = 9
        self.numberStateVariables = 30
        # This model starts with no Maxwell elements defined.
        # They are added as the material properties file is read.
        # NOTE: ** I have temporarily changed this behavior to allow exactly 3 Maxwell models.
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
