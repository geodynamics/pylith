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
      function rtsafe(x1,x2,xacc,ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,
     & efstsi,gam,dlam,deltp,alfap,iopd,idout,kto,kw,plas)
c
c...function to find the root of the effective stress function
c   adapted from Numerical Recipes
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iopd,idout,kto,kw
      double precision rtsafe,x1,x2,xacc,ae,bs,c,ds,dl1,dl2,t1,t2,emhu
      double precision anpwr,efstsi,gam,dlam,deltp,alfap
      logical plas
c
c...  defined constants
c
      include "rconsts.inc"
c
      integer maxit
      double precision eps
      parameter (maxit=100,eps=1.0d-18)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer j
      double precision xfac,acc,fl,df,fh,xl,xh,dxold,dx,f,temp
c
c*      write(6,*) "Hello from rtsafe_f!"
c
      xfac=half*(x1+x2)
      if(xfac.eq.zero) xfac=one
      acc=xacc*xfac
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x1,efstsi,
     & gam,dlam,fl,df,deltp,alfap,iopd,plas)
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x2,efstsi,
     & gam,dlam,fh,df,deltp,alfap,iopd,plas)
      if((fl.gt.zero.and.fh.gt.zero).or.(fl.lt.zero.and.fh.lt.zero))
     & pause 'Root must be bracketed in rtsafe!'
      if(fl.eq.zero) then
        rtsafe=x1
        return
      else if(fh.eq.zero) then
        rtsafe=x2
        return
      else if(fl.lt.zero)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      end if
      rtsafe=half*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,rtsafe,efstsi,
     & gam,dlam,f,df,deltp,alfap,iopd,plas)
      do 11 j=1,maxit
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.zero
     *      .or. abs(two*f).gt.abs(dxold*df) ) then
          dxold=dx
          dx=half*(xh-xl)
          rtsafe=xl+dx
          if(abs(xl-rtsafe).le.eps)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(abs(temp-rtsafe).le.eps)return
        end if
        if(abs(dx).lt.acc) return
        call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,rtsafe,efstsi,
     &   gam,dlam,f,df,deltp,alfap,iopd,plas)
        if(f.lt.zero) then
          xl=rtsafe
        else if(f.gt.zero) then
          xh=rtsafe
        else
          return
        end if
11    continue
      write(kto,*) 'rtsafe exceeding maximum iterations!'
      if(idout.gt.1) write(kw,*) 'rtsafe exceeding maximum iterations!'
      return
      end
c
c version
c $Id: rtsafe.f,v 1.1 2004/07/12 19:35:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
