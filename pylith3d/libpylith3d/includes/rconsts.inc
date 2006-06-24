c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
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
c...  rconsts_par.inc:  A set of commonly-used constants.
c
      double precision zero,one,two,three,four,five,six,seven,eight,nine
      double precision half,third,fourth,fifth,sixth,eighth,root3,root3i
      double precision pi,big
      parameter(zero   = 0.0d0,
     &          one    = 1.0d0,
     &          two    = 2.0d0,
     &          three  = 3.0d0,
     &          four   = 4.0d0,
     &          five   = 5.0d0,
     &          six    = 6.0d0,
     &          seven  = 7.0d0,
     &          eight  = 8.0d0,
     &          nine   = 9.0d0)
c
      parameter(half   = 0.5d0,
     &          third  = one/three,
     &          fourth = 0.25d0,
     &          fifth  = 0.2d0,
     &          sixth  = one/six,
     &          eighth = 0.125d0,
     &          root3  = 1.732050807568877293527446341505d0,
     &          root3i = 0.577350269189625764509148780502d0,
     &          pi     = 3.141592653589793238462643383276d0,
     &          big    = 1.0d30)
c
c version
c $Id: rconsts.inc,v 1.2 2004/05/24 20:06:21 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
