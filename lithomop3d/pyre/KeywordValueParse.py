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

class KeywordValueParse:

    def parseline(self,line):

        import pyre.units

        uparser = pyre.units.parser()

        # print "Hello from KeywordValueParse.parseline!"

        self.keyvals[2] = False
        comment = line.find('#')
        if comment == 0: return self.keyvals
        stest = line[:comment].split('=')
        if len(stest) != 2: return self.keyvals
        key = stest[0].strip()
        rawvalue = stest[1].strip()

        try:
            uvalue = uparser.parse(rawvalue)
            value =uvalue.value
        except (NameError, AttributeError):
	    try:
                value = eval(rawvalue)
	    except (NameError):
		value = rawvalue
            except:
		return self.keyvals
        except:
            return self.keyvals

        self.keyvals[0] = key
        self.keyvals[1] = value
        self.keyvals[2] = True

        return self.keyvals


    def __init__(self):
        print ""
        print "Hello from KeywordValueParse.__init__!"
        self.keyvals = [None, None, None]
        return

# version
# $Id: KeywordValueParse.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $

# End of file 
