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
      subroutine pquadwedge(sh,gauss,nen,ngauss,intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a quadratic wedge.
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
      integer it(15),ih(15)
      data it/1,2,3,1,2,3,4,5,6,4,5,6,1,2,3/
      data ih/1,1,1,2,2,2,1,1,1,2,2,2,3,3,3/
c
c...  intrinsic functions
c
      intrinsic sqrt,dble
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i,l,nshsize,ngssize
      double precision t(6),h(3),dt(6,2),dh(3)
      double precision rr,ss,tt,uu,g1,w1
      double precision trivol
c
c...  definitions
c
      trivol=third
      nshsize=(nsd+1)*nen*ngauss
      ngssize=(nsd+1)*ngauss
c
c...  quadratic wedge definition
c
      do l=1,ngauss
        gauss(1,l)=third
        gauss(2,l)=third
        gauss(3,l)=root3i
        gauss(4,l)=trivol
      end do
      gauss(3,2)=-root3i
      if(intord.ne.2) then
        do l=1,ngauss
          do i=1,nen
            g1=dble(l-2)*sqrt(three/five)
            w1=five/nine
            if(l.eq.2) w1=eight/nine
            w1=w1*trivol
            gauss(1,3*(l-1)+1)=sixth
            gauss(2,3*(l-1)+1)=sixth
            gauss(3,3*(l-1)+1)=g1
            gauss(4,3*(l-1)+1)=third*w1
            gauss(1,3*(l-1)+2)=two*third
            gauss(2,3*(l-1)+2)=sixth
            gauss(3,3*(l-1)+2)=g1
            gauss(4,3*(l-1)+2)=third*w1
            gauss(1,3*(l-1)+3)=sixth
            gauss(2,3*(l-1)+3)=two*third
            gauss(3,3*(l-1)+3)=g1
            gauss(4,3*(l-1)+3)=third*w1
          end do
        end do
      end if
c
      do l=1,ngauss
        rr=gauss(1,l)
        ss=gauss(2,l)
        tt=gauss(3,l)
        uu=one-rr-ss
        t(1)=two*rr*rr-rr
        t(2)=two*ss*ss-ss
        t(3)=two*uu*uu-uu
        t(4)=four*rr*ss
        t(5)=four*ss*uu
        t(6)=four*rr*uu
        h(1)=half*(tt*tt+tt)
        h(2)=half*(tt*tt-tt)
        h(3)=one-tt*tt
        dt(1,1)=four*rr-one
        dt(2,1)=zero
        dt(3,1)=one-four*uu
        dt(4,1)=four*ss
        dt(5,1)=(-four)*ss
        dt(6,1)=four*(uu-rr)
        dt(1,2)=zero
        dt(2,2)=four*ss-one
        dt(3,2)=one-four*uu
        dt(4,2)=four*rr
        dt(5,2)=four*(uu-ss)
        dt(6,2)=(-four)*rr
        dh(1)=half*(two*tt+one)
        dh(2)=half*(two*tt-one)
        dh(3)=one-two*tt
        do i=1,nen
          sh(4,i,l)=t(it(i))*h(ih(i))
          sh(1,i,l)=dt(it(i),1)*h(ih(i))
          sh(2,i,l)=dt(it(i),2)*h(ih(i))
          sh(3,i,l)=t(it(i))*dh(ih(i))
        end do
      end do
c
      return
      end
c
c version
c $Id: pquadwedge.f,v 1.3 2005/03/22 04:45:55 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
