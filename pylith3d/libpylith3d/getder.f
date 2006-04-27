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
      subroutine getder(det,sh,shd,xs,nen)
c
c...  subroutine to compute shape function derivatives in
c     global coordinates
c
c        sh(1,nen),sh(2,nen),sh(3,nen)    = r, s, and t derivatives
c                                           of shape functions
c        shd(1,nen),shd(2,nen),shd(3,nen) = x, y, and z derivatives
c                                           of shape functions
c        xs(nsd,nsd)                      = jacobian matrix
c        det                              = jacobian matrix determinant
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nen
      double precision det,sh(nsd+1,nen),shd(nsd+1,nen),xs(nsd,nsd)
c
c...  local variables
c
      integer i
      double precision a11,a12,a13,a21,a22,a23,a31,a32,a33,detinv
c
cdebug      write(6,*) "Hello from getder_f!"
c
c
c...transform natural derivatives to (x,y,z) derivatives using an
c   explicit 3x3 matrix inversion routine
c
cdebug      a11 = (xs(2,2)*xs(3,3))-(xs(2,3)*xs(3,2))
cdebug      a12 =-(xs(1,2)*xs(3,3))+(xs(1,3)*xs(3,2))
cdebug      a13 = (xs(1,2)*xs(2,3))-(xs(1,3)*xs(2,2))
cdebug      a21 =-(xs(2,1)*xs(3,3))+(xs(2,3)*xs(3,1))
cdebug      a22 = (xs(1,1)*xs(3,3))-(xs(1,3)*xs(3,1))
cdebug      a23 =-(xs(1,1)*xs(2,3))+(xs(1,3)*xs(2,1))
cdebug      a31 = (xs(2,1)*xs(3,2))-(xs(2,2)*xs(3,1))
cdebug      a32 =-(xs(1,1)*xs(3,2))+(xs(1,2)*xs(3,1))
cdebug      a33 = (xs(1,1)*xs(2,2))-(xs(1,2)*xs(2,1))
      a11 = (xs(2,2)*xs(3,3))-(xs(2,3)*xs(3,2))
      a12 =-(xs(2,1)*xs(3,3))+(xs(3,1)*xs(2,3))
      a13 = (xs(2,1)*xs(3,2))-(xs(3,1)*xs(2,2))
      a21 =-(xs(1,2)*xs(3,3))+(xs(3,2)*xs(1,3))
      a22 = (xs(1,1)*xs(3,3))-(xs(1,3)*xs(3,1))
      a23 =-(xs(1,1)*xs(3,2))+(xs(3,1)*xs(1,2))
      a31 = (xs(1,2)*xs(2,3))-(xs(2,2)*xs(1,3))
      a32 =-(xs(1,1)*xs(2,3))+(xs(2,1)*xs(1,3))
      a33 = (xs(1,1)*xs(2,2))-(xs(1,2)*xs(2,1))
      detinv=one/det
      do i=1,nen
        shd(1,i)=detinv*(a11*sh(1,i)+a12*sh(2,i)+a13*sh(3,i))
        shd(2,i)=detinv*(a21*sh(1,i)+a22*sh(2,i)+a23*sh(3,i))
        shd(3,i)=detinv*(a31*sh(1,i)+a32*sh(2,i)+a33*sh(3,i))
        shd(4,i)=sh(4,i)
      end do
c
      return
      end
c
c version
c $Id: getder.f,v 1.4 2004/08/12 01:25:34 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
