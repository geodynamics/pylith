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
c
      subroutine plintet(sh,gauss,nen,ngauss,intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a linear tetrahedron.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nen,ngauss,intord
      double precision sh(nsd+1,nen,ngauss)
      double precision gauss(nsd+1,ngauss)
c
c...  local constants
c
      double precision r(4),s(4),t(4),u(4)
      data u/ 1d0, 0d0, 0d0, 0d0/
      data r/ 0d0, 1d0, 0d0, 0d0/
      data s/ 0d0, 0d0, 1d0, 0d0/
      data t/ 0d0, 0d0, 0d0, 1d0/
c
c...  intrinsic functions
c
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i
      double precision rr,ss,tt,uu
      double precision tetvol
c
c...  definitions
c
      tetvol=sixth
c
c...  Linear hex definition
c     One-point integration is used in all cases.
c
      gauss(1,1)=fourth
      gauss(2,1)=fourth
      gauss(3,1)=fourth
      gauss(4,1)=tetvol
c
      do i=1,nen
        rr=r(i)*gauss(1,1)
        ss=s(i)*gauss(2,1)
        tt=t(i)*gauss(3,1)
        uu=u(i)*gauss(3,1)
        sh(4,i,1)=rr+ss+tt+uu
        sh(1,i,1)=r(i)-u(i)
        sh(2,i,1)=s(i)-u(i)
        sh(3,i,1)=t(i)-u(i)
      end do
c
      return
      end
c
c version
c $Id: plintet.f,v 1.4 2005/03/22 04:45:55 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
