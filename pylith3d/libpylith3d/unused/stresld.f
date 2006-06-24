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
      subroutine stresld(stn,scur,st0,eps,beta,dbeta,betb,dbetb,dmat,
     & ee,wr,prop,rtimdat,stol,iddmat,n,nstr,ndof,nddmat,nprop,
     & ipstrs,nprestr,nstep,lgdefp,idout,kto,kw,ivisc,iplas)
c
c...subroutine to perform stress integration for large deformations.
c   The total strain is also computed.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer n,nstr,ndof,nddmat,nprop,ipstrs,nprestr,nstep,lgdefp
      integer idout,kto,kw,ivisc,iplas
      double precision stn(nstr),scur(nstr),st0(nstr),eps(nstr)
      double precision beta(nstr),dbeta(nstr),betb(nstr),dbetb(nstr)
      double precision dmat(nddmat),ee(nstr),wr(ndof),prop(nprop)
      double precision stol
c
c...  included dimension and type statements
c
      include "iddmat_dim.inc"
      include "rtimdat_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
      integer idim
      parameter(idim=3)
c
c...  intrinsic functions
c
      intrinsic log
c
c...  local variables
c
      integer nrot,i
      double precision e1,e2,e3,rsub
      double precision r(3,3),u(3,3),ph(3,3),rtot(3,3),eet(6)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
c
c*      write(6,*) "Hello from stresld_f!"
c
      ee(4)=half*ee(4)
      ee(5)=half*ee(5)
      ee(6)=half*ee(6)
c
c...perform polar decomposition of deformation gradient
c
      call poldcmp(ee,wr,r,u)
c
c...perform eigenvalue decomposition on the stretch matrix, u, and
c   compute total rotation matrix, rtot
c
      call jacobi(u,idim,idim,ee,ph,nrot)
      if(nprestr.ne.0) then
        call dcopy(idim*idim,r,ione,rtot,ione)
      else
        call dgemm("n","n",idim,idim,idim,one,r,idim,ph,idim,zero,rtot,
     &   idim)
      end if
c
c...compute logarithmic principal stretches and logarithmic total strain
c
      e1=log(ee(1))
      e2=log(ee(2))
      e3=log(ee(3))
      eps(1)=ph(1,1)*ph(1,1)*e1+ph(1,2)*ph(1,2)*e2+ph(1,3)*ph(1,3)*e3
      eps(2)=ph(2,1)*ph(2,1)*e1+ph(2,2)*ph(2,2)*e2+ph(2,3)*ph(2,3)*e3
      eps(3)=ph(3,1)*ph(3,1)*e1+ph(3,2)*ph(3,2)*e2+ph(3,3)*ph(3,3)*e3
      eps(4)=ph(1,1)*ph(2,1)*e1+ph(1,2)*ph(2,2)*e2+ph(1,3)*ph(2,3)*e3
      eps(5)=ph(2,1)*ph(3,1)*e1+ph(2,2)*ph(3,2)*e2+ph(2,3)*ph(3,3)*e3
      eps(6)=ph(1,1)*ph(3,1)*e1+ph(1,2)*ph(3,2)*e2+ph(1,3)*ph(3,3)*e3
c*****  double-check this.  I should be computing incremental
c*****  logarithmic strains rotated to configuration at time t.
      if(((ivisc.eq.1).or.(iplas.eq.1)).and.nstep.gt.0) then
        call esfcomp(stn,scur,eps,beta,dbeta,betb,dbetb,prop,
     &   rtimdat,stol,n,nstr,ndof,nprop,ipstrs,nstep,lgdefp,idout,kto,
     &   kw,ivisc,iplas)
        eps(4)=two*eps(4)
        eps(5)=two*eps(5)
        eps(6)=two*eps(6)
      else
        eps(4)=two*eps(4)
        eps(5)=two*eps(5)
        eps(6)=two*eps(6)
        call dcopy(nstr,eps,ione,eet,ione)
        do i=1,ndof
          rsub=prop(nprop-1)
          if(rsub.gt.-1.0d0) rsub=log(one+rsub)
          eet(i)=eet(i)-third*rsub
        end do
        call dspmv("l",nstr,one,dmat,eet,ione,zero,scur,ione)
      end if
c
c...transform stresses to principal stretch directions
c
      r(1,1)=scur(1)
      r(2,2)=scur(2)
      r(3,3)=scur(3)
      r(1,2)=scur(4)
      r(2,3)=scur(5)
      r(1,3)=scur(6)
      call dsymm("l","l",idim,idim,one,r,idim,ph,idim,zero,u,idim)
      call dgemm("t","n",idim,idim,idim,one,ph,idim,u,idim,zero,r,idim)
c
c...compute Cauchy stresses in principal stretch directions
c
      u(1,1)=r(1,1)
      u(2,2)=r(2,2)
      u(3,3)=r(3,3)
      u(1,2)=r(1,2)
      u(1,3)=r(1,3)
      u(2,3)=r(2,3)
      if(ee(1).ne.ee(2)) u(1,2)=r(1,2)*(two*ee(1)*ee(2)*
     & log(ee(1)/ee(2))/(ee(1)*ee(1)-ee(2)*ee(2)))
      if(ee(1).ne.ee(3)) u(1,3)=r(1,3)*(two*ee(1)*ee(3)*
     & log(ee(1)/ee(3))/(ee(1)*ee(1)-ee(3)*ee(3)))
      if(ee(2).ne.ee(3)) u(2,3)=r(2,3)*(two*ee(2)*ee(3)*
     & log(ee(2)/ee(3))/(ee(2)*ee(2)-ee(3)*ee(3)))
      u(2,1)=u(1,2)
      u(3,1)=u(1,3)
      u(3,2)=u(2,3)
c
c...rotate to obtain final Cauchy stresses
c
      if(nprestr.ne.0) then
        call dgemm("n","n",idim,idim,idim,one,ph,idim,u,idim,zero,r,
     &   idim)
        call dgemm("n","t",idim,idim,idim,one,r,idim,ph,idim,zero,u,
     &   idim)
        u(1,1)=u(1,1)+st0(1)
        u(2,2)=u(2,2)+st0(2)
        u(3,3)=u(3,3)+st0(3)
        u(1,2)=u(1,2)+st0(4)
        u(2,3)=u(2,3)+st0(5)
        u(1,3)=u(1,3)+st0(6)
        u(2,1)=u(1,2)
        u(3,1)=u(1,3)
        u(3,2)=u(2,3)
      end if
      call dgemm("n","n",idim,idim,idim,one,rtot,idim,u,idim,zero,r,
     & idim)
      call dgemm("n","t",idim,idim,idim,one,r,idim,rtot,idim,zero,u,
     & idim)
      scur(1)=u(1,1)
      scur(2)=u(2,2)
      scur(3)=u(3,3)
      scur(4)=u(1,2)
      scur(5)=u(2,3)
      scur(6)=u(1,3)
      return
      end
c
c version
c $Id: stresld.f,v 1.1 2004/07/12 20:58:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
