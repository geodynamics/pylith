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


from pylith3d.KeywordValueParse import KeywordValueParse

class MaterialModel:


    def readprop(self,file):

        lineparse = KeywordValueParse()

        # print "Hello from MaterialModel.readprop!"

        endMaterial = False
        while not endMaterial:
            line = file.readline()
            if not line: break
            keyvals = lineparse.parseline(line)
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
        print ""
        print "Hello from MaterialModel.__init__!"
        self.propertyList = [0.0] * self.numberProperties
        return

# version
# $Id: MaterialModel.py,v 1.4 2005/01/06 02:04:21 willic3 Exp $

# End of file 
