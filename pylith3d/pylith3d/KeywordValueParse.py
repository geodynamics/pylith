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

class KeywordValueParse:

    def parseline(self,line):

        import pyre.units

        uparser = pyre.units.parser()

        self.trace.log("Hello from KeywordValueParse.parseline!")

        self.keyvals[3] = False
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
            uvalue = value
        except:
            return self.keyvals

        self.keyvals[0] = key
        self.keyvals[1] = value
        self.keyvals[2] = uvalue
        self.keyvals[3] = True

        return self.keyvals


    def __init__(self):
        import journal
        self.trace = journal.debug("pylith3d.trace")
        
        self.trace.log("Hello from KeywordValueParse.__init__!")

        self.keyvals = [None, None, None,None]
        return

# version
# $Id: KeywordValueParse.py,v 1.4 2005/01/06 01:45:13 willic3 Exp $

# End of file 
