c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c...  nshape.inc:  A header that defines element type information.
c
c     At present, the largest number of nodes per element is 20 (for a
c     quadratic hexahedron).  The largest number of Gauss points also
c     occurs for this element, with a value of 27 (for 3x3x3
c     integration).  The maximum number of element equations is derived
c     from nenmax, and the number of different element types is
c     presently 62 (this includes all of the infinite element
c     permutations).
c
c     This header should be included before any subroutine arguments,
c     as the parameters will generally be used to dimension the
c     arguments.
c
      integer nenmax,ngaussmax,neemax,netypes,netypesi
      parameter(nenmax=20,ngaussmax=27,neemax=3*nenmax,netypes=62,
     & netypesi=10)
c
c version
c $Id: nshape.inc,v 1.4 2004/07/09 16:25:17 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
