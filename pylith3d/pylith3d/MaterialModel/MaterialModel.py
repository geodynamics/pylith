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
