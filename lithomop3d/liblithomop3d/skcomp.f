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
      subroutine skcomp(x,d,skew,idslp,ipslp,nsd,ndof,nskdim,npdim,
     & ipstrs,numsn,numnp,nstep,lgdefp,kto)
c
c...  subroutine to compute skew angles for a given load vector and
c     fault orientation
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nsd,ndof,nskdim,npdim,ipstrs,numsn,numnp,nstep,lgdefp,kto
      integer idslp(numsn),ipslp(npdim,numsn)
      double precision x(nsd,numnp),d(ndof,numnp),skew(nskdim,numnp)
c
c...  defined constants
c
      include "rconsts.inc"
      double precision eps,big,rdeg
      parameter(eps=1.0d-8,big=1.0d30,rdeg=180.0d0/pi)
c
c...  intrinsic functions
c
      intrinsic acos,sign,dble,sqrt,min,sin,cos,asin,tan,atan
c
c...  local variables
c
      integer ma,nd,ncvm,ldtmp,nn,n,i
      double precision rm,sdmin,xd,yd,zd,amag,bmag,sgn,sk1
      double precision dnm,rat
      double precision a(3),xx(5),yy(5),sd(5),xl(3,5),cov(3,3)
c
c*      write(6,*) "Hello from skcomp_f!"
c
      ma=3
      nd=5
      if(numsn.eq.4) nd=4
      if(numsn.eq.3) nd=3
      if(numsn.lt.3) then
        write(kto,*) 'Not enough points to define a plane in skcomp!'
        stop
      end if
      ncvm=3
      ldtmp=lgdefp
      if(ipstrs.eq.1.and.nstep.eq.0) ldtmp=0
      rm=one
      if(ldtmp.eq.0) rm=zero
c
c...  loop over slippery nodes, and compute normal for each plane
c
      do nn=1,numsn
        n=idslp(nn)
        sdmin=big
c
c...  localize coordinates and set up values for least-squares fit to
c     a plane
c
        xl(1,nd)=x(1,n)+rm*d(1,n)
        xl(2,nd)=x(2,n)+rm*d(2,n)
        xl(3,nd)=x(3,n)+rm*d(3,n)
        xx(nd)=dble(nd)
        yy(nd)=one
        do i=1,nd-1
          xl(1,i)=x(1,ipslp(i,nn))+rm*d(1,ipslp(i,nn))
          xl(2,i)=x(2,ipslp(i,nn))+rm*d(2,ipslp(i,nn))
          xl(3,i)=x(3,ipslp(i,nn))+rm*d(3,ipslp(i,nn))
          xd=xl(1,i)-xl(1,nd)
          yd=xl(2,i)-xl(2,nd)
          zd=xl(3,i)-xl(3,nd)
          xx(i)=dble(i)
          sd(i)=sqrt(xd*xd+yd*yd+zd*zd)
          yy(i)=one
          sdmin=min(sdmin,sd(i))
        end do
        sd(nd)=half*sdmin
        call lfit(xx,yy,sd,nd,a,ma,cov,ncvm,xl)
        amag=one/sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        a(1)=a(1)*amag
        a(2)=a(2)*amag
        a(3)=a(3)*amag
        amag=one/sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        bmag=zero
        sgn=one
        if(a(3).lt.zero) sgn=-one
        do i=1,ndof
          a(i)=sgn*a(i)*amag
          if(abs(a(i)).lt.eps) a(i)=zero
          bmag=bmag+a(i)*a(i)
        end do
        bmag=one/sqrt(bmag)
        do i=1,ndof
          a(i)=a(i)*bmag
          sgn=sign(one,a(i))
          if(abs(a(i)).gt.one) a(i)=sgn
        end do
c
c...  new z-axis is perpendicular to fault plane
c
        sk1=acos(a(3))
        if(abs(a(3)-one).lt.eps) then
          skew(1,n)=zero
          skew(2,n)=zero
        else if(abs(a(1)).gt.eps) then
          skew(1,n)=atan(a(2)/a(1))
          dnm=sin(skew(1,n))+cos(skew(1,n))
          rat=(a(1)+a(2))/dnm
          sgn=sign(one,rat)
          if(abs(rat).gt.one) rat=sgn
          skew(2,n)=asin(rat)
        else if(abs(a(3)).gt.eps) then
          skew(2,n)=sk1
          dnm=a(3)*tan(skew(2,n))
          rat=a(2)/dnm
          sgn=sign(one,rat)
          if(abs(rat).gt.one) rat=sgn
          skew(1,n)=asin(rat)
        else
          skew(2,n)=half*pi
          skew(1,n)=asin(a(2))
        end if
c*        write(17,710) n,a(1),a(2),a(3),rdeg*skew(1,n),rdeg*skew(2,n)
      end do
c*      write(17,720)
c*      call flush(17)
700   format('Results for total iteration # ',i6,/)
710   format(i6,5(2x,1pe12.5))
720   format(///)
      return
      end
c
c version
c $Id: skcomp.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
