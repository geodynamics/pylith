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
      subroutine gspre(alnz,bres,z,rtz,ja,neq,nnz)

c...  subroutine to perform symmetrized Gauss-Seidel preconditioning
c     This routine returns z=B-inverse*bres, where B is an approximation
c     to the matrix alnz.  The routine also returns the dot product of
c     z and bres.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer neq,nnz
      integer ja(nnz)
      double precision alnz(nnz),bres(neq),z(neq),rtz
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  local variables
c
      integer i,j
      double precision sum,dsr
c
c...  perform forward reduction
c
      call dcopy(neq,bres,ione,z,ione)
      do i=1,neq
        sum=z(i)
        do j=ja(i),ja(i+1)-1
          if(ja(j).gt.i) go to 10
          sum=sum-z(ja(j))*alnz(j)/sqrt(alnz(ja(j)))
        end do
10      continue
        z(i)=sum/sqrt(alnz(i))
      end do
c
c...  perform back substitution
c
      rtz=zero
      do i=neq,1,-1
        dsr=one/sqrt(alnz(i))
        sum=z(i)
        do j=ja(i+1)-1,ja(i),-1
          if(ja(j).lt.i) go to 20
          sum=sum-z(ja(j))*alnz(j)*dsr
        end do
20      continue
        z(i)=sum*dsr
        rtz=rtz+bres(i)*z(i)
      end do
      return
      end
c
c version
c $Id: gspre.f,v 1.1 2005/03/30 04:51:26 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
