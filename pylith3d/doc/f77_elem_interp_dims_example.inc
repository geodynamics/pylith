c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams
c  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c...  Example include file to define element interpolation information
c     in f77.  This is a temporary template, based on a linear hex
c     element with full quadrature.  It is unclear at present whether
c     it is desirable to put all of this information into an 'include'
c     file, which will require it to be compiled as a separate element
c     type.  Since all of my present routines that deal with
c     interpolation information do not require any knowledge of the
c     element characteristics, it seems a more general approach would
c     be to either tabulate (or read from a tabulated list) this
c     information at runtime.  This would prevent the problem of
c     having to define a new element type corresponding to each possible
c     integration order, for example.
c
c...  Global dimensions, defined as parameters.
c
c       nen:     Number of element nodes
c       ngauss:  Number of element Gauss points
c
      integer nen,ngauss
      parameter(nen=8,ngauss=8)
c
c...  Note that, because of f77 requirements about statement order, it
c     is easiest if the interpolation information is broken up into 2
c     files.  This first file simply defines the global dimensions.
c
c version
c $Id: f77_elem_interp_dims_example.inc,v 1.1 2005/02/28 20:48:44 willic3 Exp $

c Generated automatically by Fortran77Mill on Mon Feb 28 11:48:22 2005

c End of file 
