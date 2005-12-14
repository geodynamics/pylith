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
      subroutine pquadpyr(sh,gauss,nen,ngauss,intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a quadratic pyramid.
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
      double precision c1,c2,c3,c4,c5
      parameter(c1 =  15.0d0,
     &          c2 =  81.0d0,
     &          c3 = 100.0d0,
     &          c4 = 125.0d0,
     &          c5 =  27.0d0)
c
      double precision g1,g2,g3,g4,g5,g6,g7,w1,w2,w3,w4
      parameter(g1 =  0.673931986207731725795838166975d0,
     &          g2 =  0.610639618865075531589990884070d0,
     &          g3 =  0.580939660561084423142479797240d0,
     &          g4 = -0.830065359477124183006535947712d0,
     &          g5 =  0.524394036075370071893249885853d0,
     &          g6 = -(one/seven),
     &          g7 = -(nine/28.0d0),
     &          w1 =  0.515003019323671497584541062801d0,
     &          w2 =  0.257183745242064658853721141309d0,
     &          w3 =  0.419515737191525949778403607347d0,
     &          w4 =  2.474004977113405935668648781299d0)
c
      double precision r(13),s(13)
      data r/-1d0, 1d0, 1d0,-1d0, 0d0, 0d0, 1d0, 0d0,-1d0,
     &       -1d0, 1d0, 1d0,-1d0/
      data s/-1d0,-1d0, 1d0, 1d0, 0d0,-1d0, 0d0, 1d0, 0d0,
     &       -1d0,-1d0, 1d0, 1d0/
c
      integer iq(13),ih(13)
      data iq/1,2,3,4,9,5,6,7,8,1,2,3,4/
      data ih/1,1,1,1,2,1,1,1,1,3,3,3,3/
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
      double precision q(9),h(3),dq(9,2),dh(3)
      double precision rr,ss,tt,g,w
c
c...  definitions
c
c
c...  Quadratic pyramid definition
c
      if(intord.eq.2) then
        g=eight*sqrt(two/c1)/five
        w=c2/c3
        do l=1,ngauss
          gauss(1,l)=r(l)*g
          gauss(2,l)=s(l)*g
          gauss(3,l)=-(two*third)
          gauss(4,l)=w
        end do
        gauss(3,5)=two/five
        gauss(4,5)=c4/c5
      else
        do l=1,4
          gauss(1,l)=r(l)*g1
          gauss(2,l)=s(l)*g1
          gauss(3,l)=g6
          gauss(4,l)=w1
        end do
        do l=5,8
          rr=zero
          ss=zero
          if(l.eq.5) rr=-one
          if(l.eq.6) rr=one
          if(l.eq.7) ss=-one
          if(l.eq.8) ss=one
          gauss(1,l)=rr*g2
          gauss(2,l)=ss*g2
          gauss(3,l)=g7
          gauss(4,l)=w2
        end do
        do l=9,12
          gauss(1,l)=r(l-8)*g3
          gauss(2,l)=s(l-8)*g3
          gauss(3,l)=g4
          gauss(4,l)=w3
        end do
        gauss(1,13)=zero
        gauss(2,13)=zero
        gauss(3,13)=g5
        gauss(4,13)=w4
      end if
c
      do l=1,ngauss
        rr=gauss(1,l)
        ss=gauss(2,l)
        tt=gauss(3,l)
        q(1)=fourth*(one-rr)*(one-ss)*(-rr-ss-one)
        q(2)=fourth*(one+rr)*(one-ss)*( rr-ss-one)
        q(3)=fourth*(one+rr)*(one+ss)*( rr+ss-one)
        q(4)=fourth*(one-rr)*(one+ss)*(-rr+ss-one)
        q(5)=half*(one-rr*rr)*(one-ss)
        q(6)=half*(one+rr)*(one-ss*ss)
        q(7)=half*(one-rr*rr)*(one+ss)
        q(8)=half*(one-rr)*(one-ss*ss)
        q(9)=one
        h(1)=half*(tt*tt-tt)
        h(2)=half*(tt*tt+tt)
        h(3)=one-tt*tt
        dq(1,1)=(-fourth)*((one-ss)*(-rr-ss-one)+(one-rr)*(one-ss))
        dq(2,1)=( fourth)*((one-ss)*( rr-ss-one)+(one+rr)*(one-ss))
        dq(3,1)=( fourth)*((one+ss)*( rr+ss-one)+(one+rr)*(one+ss))
        dq(4,1)=(-fourth)*((one+ss)*(-rr+ss-one)+(one-rr)*(one+ss))
        dq(5,1)=(-rr)*(one-ss)
        dq(6,1)=( half)*(one-ss*ss)
        dq(7,1)=(-rr)*(one+ss)
        dq(8,1)=(-half)*(one-ss*ss)
        dq(1,2)=(-fourth)*((one-rr)*(-rr-ss-one)+(one-rr)*(one-ss))
        dq(2,2)=(-fourth)*((one+rr)*( rr-ss-one)+(one+rr)*(one-ss))
        dq(3,2)=( fourth)*((one+rr)*( rr+ss-one)+(one+rr)*(one+ss))
        dq(4,2)=( fourth)*((one-rr)*(-rr+ss-one)+(one-rr)*(one+ss))
        dq(5,2)=(-half)*(one-rr*rr)
        dq(6,2)=(-ss)*(one+rr)
        dq(7,2)=( half)*(one-rr*rr)
        dq(8,2)=(-ss)*(one-rr)
        dq(9,1)=zero
        dq(9,2)=zero
        dh(1)=half*(two*tt-one)
        dh(2)=half*(two*tt+one)
        dh(3)=(-two)*tt
        do i=1,nen
          sh(4,i,l)=q(iq(i))*h(ih(i))
          sh(1,i,l)=dq(iq(i),1)*h(ih(i))
          sh(2,i,l)=dq(iq(i),2)*h(ih(i))
          sh(3,i,l)=q(iq(i))*dh(ih(i))
        end do
      end do
c
      return
      end
c
c version
c $Id: pquadpyr.f,v 1.3 2005/03/22 04:45:55 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
