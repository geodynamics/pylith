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


from lithomop3d.KeywordValueParse import KeywordValueParse

class MaterialModel:


    def readprop(self,file):

        lineparse = KeywordValueParse()

        # print "Hello from MaterialModel.readprop!"

        endMaterial = False
        while not endMaterial:
            line = file.readline()
            if not line: break
            keyvals = lineparse.parseline(line)
            if keyvals[2]:
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
# $Id: MaterialModel.py,v 1.3 2004/08/12 17:11:17 willic3 Exp $

# End of file 
