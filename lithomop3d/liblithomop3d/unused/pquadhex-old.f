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
      subroutine pquadhex(sh,shj,gauss,ngauss,nen,nsd,nenmax,ngaussmax,
     & intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a quadratic (20-node)
c    hexahedron.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nsd,nenmax,ngaussmax,intord
      integer ngauss,nen
      double precision sh(nsd+1,nenmax,ngaussmax)
      double precision shj(nsd+1,nenmax,ngaussmax)
      double precision gauss(nsd+1,ngaussmax)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
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
      double precision gquad,dgquad
c
c...  local variables
c
      integer i,l,l1,l2,l3,nshsize,ngssize
      double precision g1,w1,w2,rr,ss,tt,rrw,ssw,ttw,drr,dss,dtt,beta
      double precision betai
c
c...  function definitions
c
      gquad(beta,betai)=(half*(one-abs(betai))+half)*(one+beta*betai+
     & (abs(betai)-one)*beta*beta)
      dgquad(beta,betai)=half*betai-two*beta*(one-abs(betai))
c
c...  definitions
c
      nshsize=(nsd+1)*nenmax*ngaussmax
      ngssize=(nsd+1)*ngaussmax
c
c...  Quadratic hex definition
c
      nen=20
      ngauss=ieight
      do l=1,ngauss
        gauss(1,l)=r(l)*root3i
        gauss(2,l)=s(l)*root3i
        gauss(3,l)=t(l)*root3i
        gauss(4,l)=one
      end do
      if(intord.ne.2) then
        ngauss=27
        g1=sqrt(three/five)
        w1=five/nine
        w2=eight/nine
        l=0
        do l3=1,3
          tt=dble(l3-2)
          ttw=tt*w1+(abs(tt)-one)*w2
          do l2=1,3
            ss=dble(l2-2)
            ssw=ss*w1+(abs(ss)-one)*w2
            do l1=1,3
              rr=dble(l1-2)
              rrw=rr*w1+(abs(rr)-one)*w2
              l=l+1
              gauss(1,l)=rr*g1
              gauss(2,l)=ss*g1
              gauss(3,l)=tt*g1
              gauss(4,l)=rrw*ssw*ttw
            end do
          end do
        end do
      end if
      do l=1,ngauss
        do i=1,nen
          rr=gquad(gauss(1,l),r(i))
          ss=gquad(gauss(2,l),s(i))
          tt=gquad(gauss(3,l),t(i))
          drr=dgquad(gauss(1,l),r(i))
          dss=dgquad(gauss(2,l),s(i))
          dtt=dgquad(gauss(3,l),t(i))
          sh(4,i,l)=rr*ss*tt
          sh(1,i,l)=drr*ss*tt
          sh(2,i,l)=dss*rr*tt
          sh(3,i,l)=dtt*rr*ss
        end do
        do i=1,4
          sh(i,1,l)=sh(i,1,l)-half*(sh(i,9,l)+sh(i,12,l)+ sh(i,17,l))
          sh(i,2,l)=sh(i,2,l)-half*(sh(i,9,l)+sh(i,10,l)+sh(i,18,l))
          sh(i,3,l)=sh(i,3,l)-half*(sh(i,10,l)+sh(i,11,l)+sh(i,19,l))
          sh(i,4,l)=sh(i,4,l)-half*(sh(i,11,l)+sh(i,12,l)+sh(i,20,l))
          sh(i,5,l)=sh(i,5,l)-half*(sh(i,13,l)+sh(i,16,l)+sh(i,17,l))
          sh(i,6,l)=sh(i,6,l)-half*(sh(i,13,l)+sh(i,14,l)+sh(i,18,l))
          sh(i,7,l)=sh(i,7,l)-half*(sh(i,14,l)+sh(i,15,l)+sh(i,19,l))
          sh(i,8,l)=sh(i,8,l)-half*(sh(i,15,l)+sh(i,16,l)+sh(i,20,l))
        end do
      end do
      call dcopy(nshsize,sh,ione,shj,ione)
c
      return
      end
c
c version
c $Id: pquadhex-old.f,v 1.1 2004/07/06 19:20:28 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
