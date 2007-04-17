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


#cimport libpylith3d


def f77open(unit, file, status=None, access=None, form=None, recl=None):

    # status
    if status is None:
        status = 'unknown'
    else:
        status = status.lower()
    statusChoices = ['old', 'new', 'scratch', 'unknown']
    if not status in statusChoices:
        raise ValueError("""status='%s': 'status' must be in %s""" % (status, statusChoices))

    # access
    if access is None:
        access = 'sequential'
    else:
        access = access.lower()
    if access == 'sequential':
        # recl
        if recl is not None:
            raise ValueError("""'recl' must be omitted for 'sequential' access""")
        else:
            recl = 0
    elif access == 'direct':
        # recl
        if recl is None:
            raise ValueError("""'recl' must be given for 'direct' access""")
    else:
        raise ValueError("""access='%s': 'access' must be 'sequential' or 'direct'""" % access)
    
    if form is None:
        if access == 'sequential':
            form = 'formatted'
        else:
            form = 'unformatted'
    else:
        if not form in ['formatted', 'unformatted']:
            raise ValueError("""form='%s': 'form' must be 'formatted' or 'unformatted'""" % form)

    cdef int u, ios, rl
    u = unit
    ios = 0
    rl = recl
    libpylith3d.f77open(&u, &ios, file, status, access, form, &rl,
                        len(file), len(status), len(access), len(form))
    if ios != 0:
        raise IOError("open: Fortran I/O error %d" % ios)

    return


def f77close(unit):
    cdef int u, ios
    u = unit
    ios = 0
    libpylith3d.f77close(&u, &ios)
    if ios != 0:
        raise IOError("close: Fortran I/O error %d" % ios)
    return


# end of file 
