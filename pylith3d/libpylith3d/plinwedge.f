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
c
      subroutine plinwedge(sh,gauss,nen,ngauss,intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a linear wedge.
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
      double precision r(6),s(6),t(6),u(6)
      data r/ 1d0, 0d0, 0d0, 1d0, 0d0, 0d0/
      data s/ 0d0, 1d0, 0d0, 0d0, 1d0, 0d0/
      data t/ 1d0, 1d0, 1d0,-1d0,-1d0,-1d0/
      data u/ 0d0, 0d0, 1d0, 0d0, 0d0, 1d0/
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i,l
      double precision rr,ss,tt,uu,drr,dss,dtt
      double precision trivol
c
c...  definitions
c
      trivol=half
c
c...  Linear wedge definition
c
      gauss(1,1)=zero
      gauss(2,1)=zero
      gauss(3,1)=zero
      gauss(4,1)=two*trivol
      if(intord.ne.2) then
        do l=1,ngauss
          gauss(1,l)=third
          gauss(2,l)=third
          gauss(3,l)=root3i
          gauss(4,l)=trivol
        end do
        gauss(3,2)=-root3i
      end if 
c
      do l=1,ngauss
        do i=1,nen
          rr=one-r(i)+r(i)*gauss(1,l)
          ss=one-s(i)+s(i)*gauss(2,l)
          tt=one+t(i)*gauss(3,l)
          uu=one-u(i)+u(i)*gauss(1,l)
          drr=r(i)-u(i)
          dss=s(i)-u(i)
          dtt=t(i)
          sh(4,i,l)=half*rr*ss*tt*uu
          sh(1,i,l)=half*drr*ss*tt
          sh(2,i,l)=half*dss*rr*tt
          sh(3,i,l)=half*dtt*rr*ss
        end do
      end do
c
      return
      end
c
c version
c $Id: plinwedge.f,v 1.3 2005/03/22 04:45:55 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
