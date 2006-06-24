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
      subroutine lfit(x,y,sig,ndat,a,ma,covar,npc,xn,ierr,errstrng)
c
c...  routine to perform a weighted least-squares fit
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "rconsts.inc"
      integer mmax
      parameter (mmax=10)
c
c...  subroutine arguments
c
      integer ndat,ma,npc,ierr
      double precision x(ndat),y(ndat),sig(ndat),a(ma),covar(npc,npc)
      double precision xn(3,5)
      character errstrng*(*)
c
c...  local variables
c
      integer i,j,k
      double precision sig2i,wt,ym,beta(mmax),afunc(mmax)
c
c
cdebug      write(6,*) "Hello from lfit_f!"
c
      do j=1,ma
        do k=1,ma
          covar(j,k)=zero
        end do
        beta(j)=zero
      end do
      do i=1,ndat
        call funcs(x(i),xn,afunc,ma,ndat)
        ym=y(i)
        sig2i=one/(sig(i)*sig(i))
        do j=1,ma
          wt=afunc(j)*sig2i
          do k=1,j
            covar(k,j)=covar(k,j)+wt*afunc(k)
          end do
          beta(j)=beta(j)+ym*wt
        end do
      end do
      call choldc2(covar,npc,ma,a,ierr,errstrng)
      call cholsl(covar,npc,ma,a,beta,beta)
      do j=1,ma
        a(j)=beta(j)
      end do
      return
      end
c
c version
c $Id: lfit.f,v 1.3 2004/08/12 01:39:52 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
