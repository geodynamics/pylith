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
      subroutine pquadhex(sh,gauss,nen,ngauss,intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a quadratic (20-node)
c    hexahedron.
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
      double precision r(20),s(20),t(20)
      data r/-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0,
     &        0d0, 1d0, 0d0,-1d0, 0d0, 1d0, 0d0,-1d0,
     &       -1d0, 1d0, 1d0,-1d0/
      data s/-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,
     &       -1d0, 0d0, 1d0, 0d0,-1d0, 0d0, 1d0, 0d0,
     &       -1d0,-1d0, 1d0, 1d0/
      data t/ 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0,
     &        1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0,
     &        0d0, 0d0, 0d0, 0d0/
c
c...  intrinsic functions
c
      intrinsic abs,sqrt,dble
c
c...  user-defined functions
c
      double precision gquad,dgquad,hquad,dhquadr,dhquads,dhquadt
c
c...  local variables
c
      integer i,l,l1,l2,l3,nshsize,ngssize
      double precision g1,rr,ss,tt,rrw,ssw,ttw,drr,dss,dtt
      double precision uu,dur,dus,dut,v,vi,rp,ri,sp,si,tp,ti
      double precision w(3)
cdebug      integer idb,jdb,kdb
c
c...  function definitions
c
      gquad(v,vi)=(half*(one-abs(vi))+half)*(one+v*vi+(abs(vi)-one)*v*v)
      dgquad(v,vi)=(half*(one-abs(vi))+half)*(vi+two*v*(abs(vi)-one))
      hquad(rp,sp,tp,ri,si,ti)=abs(ri*si*ti)*(ri*rp+si*sp+ti*tp-two)+
     & one-abs(ri*si*ti)
      dhquadr(ri,si,ti)=abs(ri*si*ti)*ri
      dhquads(ri,si,ti)=abs(ri*si*ti)*si
      dhquadt(ri,si,ti)=abs(ri*si*ti)*ti
c
cdebug      write(6,*) "Hello from pquadhex_f!"
c
c
c...  definitions
c
      nshsize=(nsd+1)*nen*ngauss
      ngssize=(nsd+1)*ngauss
c
c...  Quadratic hex definition
c
      do l=1,ngauss
        gauss(1,l)=r(l)*root3i
        gauss(2,l)=s(l)*root3i
        gauss(3,l)=t(l)*root3i
        gauss(4,l)=one
      end do
      if(intord.ne.2) then
        g1=sqrt(three/five)
        w(1)=five/nine
        w(2)=eight/nine
        w(3)=five/nine
        l=0
        do l3=1,3
          tt=dble(l3-itwo)
          ttw=w(l3)
          do l2=1,3
            ss=dble(l2-itwo)
            ssw=w(l2)
            do l1=1,3
              rr=dble(l1-itwo)
              rrw=w(l1)
              l=l+1
              gauss(1,l)=rr*g1
              gauss(2,l)=ss*g1
              gauss(3,l)=tt*g1
              gauss(4,l)=rrw*ssw*ttw
            end do
          end do
        end do
      end if
c
      do l=1,ngauss
        do i=1,nen
          rr=gquad(gauss(1,l),r(i))
          ss=gquad(gauss(2,l),s(i))
          tt=gquad(gauss(3,l),t(i))
          uu=hquad(gauss(1,l),gauss(2,l),gauss(3,l),r(i),s(i),t(i))
          drr=dgquad(gauss(1,l),r(i))
          dss=dgquad(gauss(2,l),s(i))
          dtt=dgquad(gauss(3,l),t(i))
          dur=dhquadr(r(i),s(i),t(i))
          dus=dhquads(r(i),s(i),t(i))
          dut=dhquadt(r(i),s(i),t(i))
          sh(4,i,l)=rr*ss*tt*uu
          sh(1,i,l)=drr*ss*tt*uu+dur*rr*ss*tt
          sh(2,i,l)=dss*rr*tt*uu+dus*rr*ss*tt
          sh(3,i,l)=dtt*rr*ss*uu+dut*rr*ss*tt
        end do
      end do
c
      return
      end
c
c version
c $Id: pquadhex.f,v 1.6 2005/03/22 04:45:55 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
