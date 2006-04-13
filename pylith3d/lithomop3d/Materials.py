#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
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
# $Id: Materials.py,v 1.3 2004/08/12 17:11:17 willic3 Exp $

# End of file 
