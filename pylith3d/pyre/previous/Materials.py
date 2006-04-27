#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams
#  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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

from MaterialModel import *

class Materials:

    def readprop(self, propertyFile):

        # print "Hello from Materials.readprop!"

        file = open(propertyFile, 'r')
        while 1:
            line = file.readline()
	    # print line
            if not line: break
            materialType = None
            exec line
            if materialType != None:
                matchoose = 'matmodel = ' + materialType + '.' + materialType + '()'
                exec matchoose
                matmodel.readprop(file)
                self.numberMaterials += 1
                self.materialNumber += [self.numberMaterials]
                self.materialModel += [matmodel.materialModel]
                self.numberProperties += [matmodel.numberProperties]
                self.propertyIndex += [len(self.propertyList) + 1]
                self.propertyList += matmodel.propertyList
		# print self.numberMaterials
		# print self.materialNumber
		# print self.materialModel
		# print self.numberProperties
		# print self.propertyIndex
		# print self.propertyList
        return self.numberMaterials


    def __init__(self):
	print ""
        print "Hello from Materials.__init__!"
        self.numberMaterials = 0
        self.materialNumber = []
        self.materialModel = []
        self.numberProperties = []
        self.propertyList = []
        self.propertyIndex = []
        return

# version
# $Id: Materials.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $

# End of file 
