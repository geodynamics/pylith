c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
c  All worldwide rights reserved.  A license to use, copy, modify and
c  distribute this software for non-commercial research purposes only
c  is hereby granted, provided that this copyright notice and
c  accompanying disclaimer is not modified or removed from the software.
c
c  DISCLAIMER:  The software is distributed "AS IS" without any express
c  or implied warranty, including but not limited to, any implied
c  warranties of merchantability or fitness for a particular purpose
c  or any warranty of non-infringement of any current or pending patent
c  rights.  The authors of the software make no representations about
c  the suitability of this software for any particular purpose.  The
c  entire risk as to the quality and performance of the software is with
c  the user.  Should the software prove defective, the user assumes the
c  cost of all necessary servicing, repair or correction.  In
c  particular, neither Rensselaer Polytechnic Institute, nor the authors
c  of the software are liable for any indirect, special, consequential,
c  or incidental damages related to the software, to the maximum extent
c  the law permits.
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
