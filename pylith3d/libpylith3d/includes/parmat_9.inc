c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2007 Rensselaer Polytechnic Institute
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
c...  parmat_9.inc:  Parameter definitions for material type 9.
c
c...  The rpar array contains:
c     1-6:    Constant part of vector, which includes current total
c             strain, viscous strain from previous step, stress-related
c             quantities from previous time step, and prestress term.
c     7-27:   Inverted elastic material matrix.
c     28:     Integration parameter alfap.
c     29:     Time step size deltp.
c     30:     Viscosity coefficient emhu.
c     31:     Power-law exponent anpwr.
c
c...  The ipar array contains:
c     1-36:   The iddmat index array.
c
      integer indconst,inddmati,indalfap,inddeltp,indemhu,indanpwr
      integer indiddmat
c
      parameter(indconst=1,
     &          inddmati=7,
     &          indalfap=28,
     &          inddeltp=29,
     &          indemhu=30,
     &          indanpwr=31,
     &          indiddmat=1)
c
c version
c $Id: parmat_9.inc,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
