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
      subroutine lfit(x,y,sig,ndat,a,ma,covar,npc,xn)
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
      integer ndat,ma,npc
      double precision x(ndat),y(ndat),sig(ndat),a(ma),covar(npc,npc)
      double precision xn(3,5)
c
c...  local variables
c
      integer i,j,k
      double precision sig2i,wt,ym,beta(mmax),afunc(mmax)
c
c
c*      write(6,*) "Hello from lfit_f!"
c
      do j=1,ma
        do k=1,ma
          covar(j,k)=zero
        end do
        beta(j)=zero
      end do
      do i=1,ndat
        call funcs(x(i),xn,afunc,ma)
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
      call choldc2(covar,npc,ma,a)
      call cholsl(covar,npc,ma,a,beta,beta)
      do j=1,ma
        a(j)=beta(j)
      end do
      return
      end
c
c version
c $Id: lfit.f,v 1.2 2004/07/07 15:55:35 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
