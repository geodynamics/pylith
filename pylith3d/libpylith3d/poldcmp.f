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
      subroutine poldcmp(ee,wr,r,u)
c
c...routine to perform a polar decomposition on the deformation
c   gradient, given the linear strain and rotation matrices
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "rconsts.inc"
c
      integer idim
      parameter(idim=3)
c
c...  subroutine arguments
c
      double precision ee(6),wr(3),r(3,3),u(3,3)
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  local variables
c
      integer nrot
      double precision e1,e2,e3
      double precision c(3,3),x(3,3)
c
c...form inverse deformation gradient and inverse Cauchy deformation
c   tensor
c
cdebug      write(6,*) "Hello from poldcmp_f!"
c
      x(1,1)=one-ee(1)
      x(2,2)=one-ee(2)
      x(3,3)=one-ee(3)
      x(1,2)=-ee(4)+wr(1)
      x(2,1)=-ee(4)-wr(1)
      x(2,3)=-ee(5)+wr(3)
      x(3,2)=-ee(5)-wr(3)
      x(1,3)=-ee(6)+wr(2)
      x(3,1)=-ee(6)-wr(2)
      call dsyrk("u","n",idim,idim,one,x,idim,zero,c,idim)
      call symmet(c,idim)
c
c...perform eigenvalue decomposition
c
c******  need to use a LAPACK eigenvalue routine -- then symmetrization
c******  probably won't be necessary.
      call jacobi(c,idim,idim,ee,r,nrot)
c
c...compute stretch tensor and rotation matrix
c
      e1=zero
      e2=zero
      e3=zero
      if(ee(1).ne.zero) e1=one/sqrt(ee(1))
      if(ee(2).ne.zero) e2=one/sqrt(ee(2))
      if(ee(3).ne.zero) e3=one/sqrt(ee(3))
      u(1,1)=r(1,1)*r(1,1)*e1+r(1,2)*r(1,2)*e2+r(1,3)*r(1,3)*e3
      u(2,2)=r(2,1)*r(2,1)*e1+r(2,2)*r(2,2)*e2+r(2,3)*r(2,3)*e3
      u(3,3)=r(3,1)*r(3,1)*e1+r(3,2)*r(3,2)*e2+r(3,3)*r(3,3)*e3
      u(1,2)=r(1,1)*r(2,1)*e1+r(1,2)*r(2,2)*e2+r(1,3)*r(2,3)*e3
      u(1,3)=r(1,1)*r(3,1)*e1+r(1,2)*r(3,2)*e2+r(1,3)*r(3,3)*e3
      u(2,3)=r(2,1)*r(3,1)*e1+r(2,2)*r(3,2)*e2+r(2,3)*r(3,3)*e3
      u(2,1)=u(1,2)
      u(3,1)=u(1,3)
      u(3,2)=u(2,3)
      call dgemm("t","t",idim,idim,idim,one,x,idim,u,idim,zero,r,idim)
      return
      end
c
c version
c $Id: poldcmp.f,v 1.4 2005/04/08 00:41:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
