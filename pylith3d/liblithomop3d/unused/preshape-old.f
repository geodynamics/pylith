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
      subroutine preshape(sh,shj,gauss,ngauss,nen,nsd,nenmax,ngaussmax,
     & netypes,intord)
c
c...program to compute the shape functions in natural coordinates.
c   The program also assigns integration points and weights.
c   All shape functions are stored in a single array.
c
c   The arguments are:
c     Input:
c
c       nsd                                = number of spatial
c                                            dimensions
c       nenmax                             = maximum number of nodes
c                                            per element
c       ngaussmax                          = maximum number of gauss
c                                            points per element
c       netypes                            = number of element types
c                                            (currently equal to 62)
c       intord                             = integration order
c                                            1 = full
c                                            2 = reduced
c                                            3 = Bbar
c
c     Output:
c       sh(nsd+1,nenmax,
c          ngaussmax,netypes)   = shape functions and their derivatives
c       shj(nsd+1,nenmax,
c           ngaussmax,netypes)  = shape functions and derivatives used
c                                 for calculating Jacobian.  Unless the
c                                 element is an infinite element, this
c                                 is identical to sh.
c       gauss(nsd+1,ngaussmax,netypes)     = Gauss point coordinates and
c                                            weights
c       ngauss(netypes)                    = number of Gauss points
c       nen(netypes)                       = number of element nodes
c
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nsd,nenmax,ngaussmax,netypes,intord
      integer ngauss(netypes),nen(netypes)
      real*8 sh(nsd+1,nenmax,ngaussmax,netypes)
      real*8 shj(nsd+1,nenmax,ngaussmax,netypes)
      real*8 gauss(nsd+1,ngaussmax,netypes)
c
c...  data statement definitions
c
      integer izero,ione,itwo,ithree,ifour,ifive,isix,iseven,ieight
      integer inine
      data izero,ione,itwo,ithree,ifour,ifive,isix,iseven,ieight,inine/
     & 0,1,2,3,4,5,6,7,8,9/
c
      real*8 zero,one,two,three,four,five,six,seven,eight,nine
      data zero,one,two,three,four,five,six,seven,eight,nine/
     & 0d0,1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0/
c
      integer idiv(3)
      data idiv/1,10,100/
c
      real*8 rlh(8),slh(8),tlh(8)
      data rlh/-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0/
      data slh/-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0/
      data tlh/ 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0/
c
      real*8 rlw(6),slw(6),tlw(6),ulw(6)
      data rlw/ 1d0, 0d0,-1d0, 1d0, 0d0,-1d0/
      data slw/ 0d0, 1d0,-1d0, 0d0, 1d0,-1d0/
      data tlw/ 1d0, 1d0, 1d0,-1d0,-1d0,-1d0/
      data ulw/ 0d0, 0d0, 1d0, 0d0, 0d0, 1d0/
c
      real*8 rlp(5),slp(5),tlp(5)
      data rlp/-1d0, 1d0, 1d0,-1d0, 0d0/
      data slp/-1d0,-1d0, 1d0, 1d0, 0d0/
      data tlp/-1d0,-1d0,-1d0,-1d0, 1d0/
c
      real*8 rlt(4),slt(4),tlt(4),ult(4)
      data rlt/ 1d0, 0d0, 0d0,-1d0/
      data slt/ 0d0, 1d0, 0d0,-1d0/
      data tlt/ 0d0, 0d0, 1d0,-1d0/
      data ult/ 0d0, 0d0, 0d0, 1d0/
c
      real*8 rqh(20),sqh(20),tqh(20)
      data rqh/-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0,
     &          0d0, 1d0, 0d0,-1d0, 0d0, 1d0, 0d0,-1d0,
     &         -1d0, 1d0, 1d0,-1d0/
      data sqh/-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,
     &         -1d0, 0d0, 1d0, 0d0,-1d0, 0d0, 1d0, 0d0,
     &         -1d0,-1d0, 1d0, 1d0/
      data tqh/ 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0,
     &          1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0,
     &          0d0, 0d0, 0d0, 0d0/
c
      real*8 rqp(13),sqp(13),tqp(13)
      data rqp/-1d0, 1d0, 1d0,-1d0, 0d0, 0d0, 1d0, 0d0,-1d0,
     &         -1d0, 1d0, 1d0,-1d0/
      data sqp/-1d0,-1d0, 1d0, 1d0, 0d0,-1d0, 0d0, 1d0, 0d0,
     &         -1d0,-1d0, 1d0, 1d0/
      data tqp/-1d0,-1d0,-1d0,-1d0, 1d0,-1d0,-1d0,-1d0,-1d0,
     &          0d0, 0d0, 0d0, 0d0/
c
      integer iqwt(15),iqwh(15)
      data iqwt/1,2,3,1,2,3,4,5,6,4,5,6,1,2,3/
      data iqwh/1,1,1,2,2,2,1,1,1,2,2,2,3,3,3/
      real*8 t(6),h(3),dt(6,2),dh(3),q(8),dq(8,2)
c
      integer iqpq(13),iqph(13)
      data iqpq/1,2,3,4,9,5,6,7,8,1,2,3,4/
      data iqph/1,1,1,1,2,1,1,1,1,3,3,3,3/
      real*8 t(6),h(3),dt(6,2),dh(3),q(9),dq(9,2)
c
c...  constants
c
      real*8 half,third,fourth,sixth,eighth,root3,root3i,trivol,tetvol
      real*8 g1,g2,g3,g4,g5,w1,w2,w3,w4,w5
      real*8 c1,c2,c3,c4,c5
      integer i,l,n,nshsize,ngssize
      real*8 rr,ss,tt,uu,drr,dss,dtt,duu
      real*8 rrw,ssw,ttw
      real*8 gquad,dgquad,beta,betai
      real*8 a11,a12,a13,a21,a22,a23,a31,a32,a33,detinv
      real*8 xs(3,3),sht(3),ra(8),sa(8),ta(8),shi(4,8)
      integer idiv(3),io(3)
c
c...  functions to compute shape functions for a quadratic hexahedron
c
      gquad(beta,betai)=(half*(one-abs(betai))+half)*(one+beta*betai+
     & (abs(betai)-one)*beta*beta)
      dgquad(beta,betai)=half*betai-two*beta*(one-abs(betai))
c
c...  definitions
c
      half=0.5d0
      third=one/three
      sixth=half*third
      fourth=0.25d0
      eighth=0.125d0
      root3=sqrt(three)
      root3i=one/root3
      trivol=half
      tetvol=sixth
      nshsize=(nsd+1)*nenmax*ngaussmax
      ngssize=(nsd+1)*ngaussmax
      call iclear(nen,netypes)
      call iclear(ngauss,netypes)
      call clear(gauss,ngssize*netypes)
      call clear(sh,nshsize*netypes)
c
c... First type:  linear hex
c
      nen(1)=ieight
      ngauss(1)=ione
      gauss(1,1,1)=zero
      gauss(2,1,1)=zero
      gauss(3,1,1)=zero
      gauss(4,1,1)=eight
      if(intord.ne.2) then
        ngauss(1)=ieight
        do l=1,ngauss(1)
          gauss(1,l,1)=rlh(l)*root3i
          gauss(2,l,1)=slh(l)*root3i
          gauss(3,l,1)=tlh(l)*root3i
          gauss(4,l,1)=one
        end do
      end if
      do l=1,ngauss(1)
        do i=1,nen(1)
          rr=one+rlh(i)*gauss(1,l,1)
          ss=one+slh(i)*gauss(2,l,1)
          tt=one+tlh(i)*gauss(3,l,1)
          sh(4,i,l,1)=eighth*rr*ss*tt
          sh(1,i,l,1)=eighth*ra(i)*ss*tt
          sh(2,i,l,1)=eighth*sa(i)*rr*tt
          sh(3,i,l,1)=eighth*ta(i)*rr*ss
        end do
      end do
      call dcopy(nshsize,sh,ione,shj,ione)
c
c...  Types 2-27:  linear hex + infinite boundaries
c     Use mod3 arithmetic to cover all cases
c
      do n=2,27
        nen(n)=nen(1)
        ngauss(n)=ngauss(1)
        call dcopy(ngssize,gauss(1,1,1),ione,gauss(1,1,n),ione)
        call dcopy(nshsize,sh(1,1,1,1),ione,sh(1,1,1,n),ione)
        io(3)=(n-1)/inine
        io(2)=(n-1-io(3)*inine)/ithree
        io(1)=mod(n-1,ithree)
        do l=1,ngauss(n)
          do i=1,nen(n)
            call infellh(gauss(1,l,n),rlh(i),rr,drr,io(1))
            call infellh(gauss(2,l,n),slh(i),ss,dss,io(2))
            call infellh(gauss(3,l,n),tlh(i),tt,dtt,io(3))
            shj(4,i,l,n)=rr*ss*tt
            shj(1,i,l,n)=drr*ss*tt
            shj(2,i,l,n)=dss*rr*tt
            shj(3,i,l,n)=dtt*rr*ss
          end do
        end do
      end do
c
c...  Type 28:  linear hex with one set of collapsed nodes (wrick)
c
      nen(28)=7
      ngauss(28)=ngauss(1)
      call dcopy(ngssize,gauss(1,1,1),ione,gauss(1,1,28),ione)
      call dcopy(nshsize,sh(1,1,1,1),ione,sh(1,1,1,28),ione)
      do l=1,ngauss(28)
        do i=1,4
          sh(i,7,l,28)=sh(i,7,l,28)+sh(i,8,l,28)
          sh(i,8,l,28)=zero
        end do
      end do
      call dcopy(nshsize,sh(1,1,1,28),ione,shj(1,1,1,28),ione)
c
c...  Type 29:  linear hex with two sets of collapsed nodes (wedge)
c     r and s are triangular coordinates in cross-section, and t is
c     cartesian natural coordinate along the length of the wedge.
c
      nen(29)=6
      ngauss(29)=ione
      gauss(1,1,29)=third
      gauss(2,1,29)=third
      gauss(3,1,29)=zero
      gauss(4,1,29)=two*trivol
      if(intord.ne.2) then
        ngauss(29)=itwo
        gauss(1,1,29)=third
        gauss(1,2,29)=third
        gauss(2,1,29)=third
        gauss(2,2,29)=third
        gauss(3,1,29)=root3i
        gauss(3,2,29)=-root3i
        gauss(4,1,29)=trivol
        gauss(4,2,29)=trivol
      end if
      do l=1,ngauss(29)
        do i=1,nen(29)
          rr=rlw(i)*gauss(1,l,29)
          ss=slw(i)*gauss(2,l,29)
          tt=one+tlw(i)*gauss(3,l,29)
          uu=half*(ulw(i)+rr+ss)
          sh(4,i,l,29)=uu*tt
          sh(1,i,l,29)=half*rlw(i)*tt
          sh(2,i,l,29)=half*slw(i)*tt
          sh(3,i,l,29)=tlw(i)*uu
        end do
      end do
      call dcopy(nshsize,sh(1,1,1,29),ione,shj(1,1,1,29),ione)
c
c...  Type 30:  linear hex with 4 points collapsed to a point (pyramid)
c     r and s are cartesian natural coordinates on the base, and t is
c     a coordinate going from negative one at the center of the base to
c     one at the apex.
c
      nen(30)=5
      ngauss(30)=ione
      gauss(1,1,30)=zero
      gauss(2,1,30)=zero
      gauss(3,1,30)=-half
      c1=128.0d0
      c2=27.0d0
      gauss(4,1,30)=c1/c2
      if(intord.ne.2) then
        c3=15.0d0
        g1=eight*sqrt(two/c3)/five
        c4=81.0d0
        c5=100.0d0
        w1=c4/c5
        ngauss(30)=ifive
        do l=1,ngauss(30)
          gauss(1,l,30)=rlp(l)*g1
          gauss(2,l,30)=slp(l)*g1
          gauss(3,l,30)=-two*third
          gauss(4,l,30)=w1
        end do
        gauss(3,5,30)=two/five
        c1=125.0d0
        gauss(4,5,30)=c1/c2
      end if
      do l=1,ngauss(30)
        do i=1,nen(30)
          rr=half*(one+rlp(i)*gauss(1,l,30))
          ss=half*(one+slp(i)*gauss(2,l,30))
          tt=half*(one-tlp(i)*gauss(3,l,30))
          sh(4,i,l,30)=rr*ss*tt
          sh(1,i,l,30)=rlp(i)*ss*tt
          sh(2,i,l,30)=slp(i)*rr*tt
          sh(3,i,l,30)=tlp(i)*rr*ss
        end do
        sh(3,nen(30),l,30)=two*sh(3,nen(30),l,30)
      end do
      call dcopy(nshsize,sh(1,1,1,30),ione,shj(1,1,1,30),ione)
c
c...  Type 31:  linear tetrahedron
c     r, s, and t are tetrahedral coordinates.
c     One-point integration is used in all cases.
c
      nen(31)=4
      ngauss(31)=ione
      gauss(1,1,31)=fourth
      gauss(2,1,31)=fourth
      gauss(3,1,31)=fourth
      gauss(4,1,31)=tetvol
      do l=1,ngauss(31)
        do i=1,nen(31)
          rr=rlt(i)*gauss(1,l,31)
          ss=slt(i)*gauss(2,l,31)
          tt=tlt(i)*gauss(3,l,31)
          uu=ult(i)+rr+ss+tt
          sh(4,i,l,31)=uu
          sh(1,i,l,31)=rlt(i)
          sh(2,i,l,31)=slt(i)
          sh(3,i,l,31)=tlt(i)
        end do
      end do
      call dcopy(nshsize,sh(1,1,1,31),ione,shj(1,1,1,31),ione)
c
c... Type 32:  quadratic (20-node) hex
c
      nen(32)=20
      ngauss(32)=ieight
      do l=1,ngauss(32)
        gauss(1,l,32)=rlh(l)*root3i
        gauss(2,l,32)=slh(l)*root3i
        gauss(3,l,32)=tlh(l)*root3i
        gauss(4,l,32)=one
      end do
      if(intord.ne.2) then
        ngauss(32)=27
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
              gauss(1,l,32)=rr*g1
              gauss(2,l,32)=ss*g1
              gauss(3,l,32)=tt*g1
              gauss(4,l,32)=rrw*ssw*ttw
            end do
          end do
        end do
      end if
      do l=1,ngauss(32)
        do i=1,nen(32)
          rr=gquad(gauss(1,l,32),rqh)
          ss=gquad(gauss(2,l,32),sqh)
          tt=gquad(gauss(3,l,32),tqh)
          drr=dgquad(gauss(1,l,32),rqh)
          dss=dgquad(gauss(2,l,32),sqh)
          dtt=dgquad(gauss(3,l,32),tqh)
          sh(4,i,l,32)=rr*ss*tt
          sh(1,i,l,32)=drr*ss*tt
          sh(2,i,l,32)=dss*rr*tt
          sh(3,i,l,32)=dtt*rr*ss
        end do
        do i=1,4
          sh(i,1,l,32)=sh(i,1,l,32)-half*(sh(i,9,l,32)+sh(i,12,l,32)+
     &     sh(i,17,l,32))
          sh(i,2,l,32)=sh(i,2,l,32)-half*(sh(i,9,l,32)+sh(i,10,l,32)+
     &     sh(i,18,l,32))
          sh(i,3,l,32)=sh(i,3,l,32)-half*(sh(i,10,l,32)+sh(i,11,l,32)+
     &     sh(i,19,l,32))
          sh(i,4,l,32)=sh(i,4,l,32)-half*(sh(i,11,l,32)+sh(i,12,l,32)+
     &     sh(i,20,l,32))
          sh(i,5,l,32)=sh(i,5,l,32)-half*(sh(i,13,l,32)+sh(i,16,l,32)+
     &     sh(i,17,l,32))
          sh(i,6,l,32)=sh(i,6,l,32)-half*(sh(i,13,l,32)+sh(i,14,l,32)+
     &     sh(i,18,l,32))
          sh(i,7,l,32)=sh(i,7,l,32)-half*(sh(i,14,l,32)+sh(i,15,l,32)+
     &     sh(i,19,l,32))
          sh(i,8,l,32)=sh(i,8,l,32)-half*(sh(i,15,l,32)+sh(i,16,l,32)+
     &     sh(i,20,l,32))
        end do
      end do
      call dcopy(nshsize,sh(1,1,1,32),ione,shj(1,1,1,32),ione)
c
c...  Types 33-58:  quadratic hex + infinite boundaries
c     Use mod3 arithmetic to cover all cases
c
      do n=33,58
        nen(n)=nen(32)
        ngauss(n)=ngauss(32)
        call dcopy(ngssize,gauss(1,1,32),ione,gauss(1,1,n),ione)
        call dcopy(nshsize,sh(1,1,1,32),ione,sh(1,1,1,n),ione)
        io(3)=(n-1)/inine
        io(2)=(n-1-io(3)*inine)/ithree
        io(1)=mod(n-1,ithree)
        do l=1,ngauss(n)
          do i=1,nen(n)
            call infelqh(gauss(1,l,n),rqh(i),rr,drr,io(1))
            call infelqh(gauss(2,l,n),sqh(i),ss,dss,io(2))
            call infelqh(gauss(3,l,n),tqh(i),tt,dtt,io(3))
            shj(4,i,l,n)=rr*ss*tt
            shj(1,i,l,n)=drr*ss*tt
            shj(2,i,l,n)=dss*rr*tt
            shj(3,i,l,n)=dtt*rr*ss
          end do
        end do
      end do
c
c...  Type 59:  quadratic hex with one set of collapsed nodes (wrick)
c
      nen(59)=18
      ngauss(59)=ngauss(32)
      call dcopy(ngssize,gauss(1,1,32),ione,gauss(1,1,59),ione)
      call dcopy(nshsize,sh(1,1,1,32),ione,sh(1,1,1,59),ione)
      do l=1,ngauss(59)
        do i=1,4
          sh(i,7,l,59)=sh(i,7,l,59)+sh(i,8,l,59)+sh(i,15,l,59)
          sh(i,8,l,59)=zero
          sh(i,15,l,59)=zero
        end do
      end do
      call dcopy(nshsize,sh(1,1,1,59),ione,shj(1,1,1,59),ione)
c
c...  Type 60:  quadratic hex with three sets of collapsed nodes (wedge)
c     r and s are triangular coordinates in cross-section, and t is
c     cartesian natural coordinate along the length of the wedge.
c
      nen(60)=15
      ngauss(60)=itwo
      gauss(1,1,60)=third
      gauss(1,2,60)=third
      gauss(2,1,60)=third
      gauss(2,2,60)=third
      gauss(3,1,60)=root3i
      gauss(3,2,60)=-root3i
      gauss(4,1,60)=trivol
      gauss(4,2,60)=trivol
      if(intord.ne.2) then
        ngauss(60)=9
        do l=1,3
          g1=dble(l-2)*sqrt(three/five)
          w1=five/nine
          if(l.eq.2) w1=eight/nine
          w1=w1*trivol
          gauss(1,3*(l-1)+1,60)=sixth
          gauss(2,3*(l-1)+1,60)=sixth
          gauss(3,3*(l-1)+1,60)=g1
          gauss(4,3*(l-1)+1,60)=third*w1
          gauss(1,3*(l-1)+2,60)=two*third
          gauss(2,3*(l-1)+2,60)=sixth
          gauss(3,3*(l-1)+2,60)=g1
          gauss(4,3*(l-1)+2,60)=third*w1
          gauss(1,3*(l-1)+3,60)=sixth
          gauss(2,3*(l-1)+3,60)=two*third
          gauss(3,3*(l-1)+3,60)=g1
          gauss(4,3*(l-1)+3,60)=third*w1
        end do
      end if
      do l=1,ngauss(60)
        rr=gauss(1,l,60)
        ss=gauss(2,l,60)
        tt=gauss(3,l,60)
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
        dt(5,1)=-four*ss
        dt(6,1)=four*(uu-rr)
        dt(1,2)=zero
        dt(2,2)=four*ss-one
        dt(3,2)=one-four*uu
        dt(4,2)=four*rr
        dt(5,2)=four*(uu-ss)
        dt(6,2)=-four*rr
        dh(1)=half*(two*tt+one)
        dh(2)=half*(two*tt-one)
        dh(3)=one-two*tt
        do i=1,nen(6)
          sh(4,i,l,60)=t(iqwt(i))*h(iqwh(i))
          sh(1,i,l,60)=dt(iqwt(i),1)*h(iqwh(i))
          sh(2,i,l,60)=dt(iqwt(i),2)*h(iqwh(i))
          sh(3,i,l,60)=t(iqwt(i))*dh(iqwh(i))
        end do
      end do
      call dcopy(nshsize,sh(1,1,1,60),ione,shj(1,1,1,60),ione)
c
c...  Type 61:  quadratic hex with 9 points collapsed to a point
c     (pyramid).  r and s are cartesian natural coordinates on the base,
c     and t is a coordinate going from negative one at the center of the
c     base to one at the apex.
c
      nen(61)=13
      c1=15.0d0
      c2=81.0d0
      c3=100.0d0
      c4=125.0d0
      c5=27.0d0
      rr=eight*sqrt(two/c1)/five
      ss=c2/c3
      ngauss(61)=ifive
      do l=1,ngauss(30)
        gauss(1,l,61)=rlp(l)*rr
        gauss(2,l,61)=slp(l)*rr
        gauss(3,l,61)=-two*third
        gauss(4,l,61)=ss
      end do
      gauss(3,5,30)=two/five
      gauss(4,5,30)=c4/c5
      if(intord.ne.2) then
        ngauss(61)=13
c*        g1=7.0*sqrt(35.0/59.0)/8.0
c*        g2=224.0*sqrt(336633710.0/33088740423.0)/37.0
c*        g3=sqrt(37043.0/35.0)/56.0
c*        g4=-127.0/153.0
c*        g5=1490761.0/2842826.0
c*        w1=170569.0/331200.0
c*        w2=276710106577408.0/1075923777052725.0
c*        w3=12827693806929.0/30577384040000.0
c*        w4=10663383340655070643544192.0/4310170528879365193704375.0
        g1=0.673931986207731725795838166975d0
        g2=0.610639618865075531589990884070d0
        g3=0.580939660561084423142479797240d0
        g4=-0.830065359477124183006535947712d0
        g5=0.524394036075370071893249885853d0
        g6=-one/seven
        g7=-nine/28.0d0
        w1=0.515003019323671497584541062801d0
        w2=0.257183745242064658853721141309d0
        w3=0.419515737191525949778403607347d0
        w4=2.474004977113405935668648781299d0
        do l=1,4
          gauss(1,l,61)=rqp(l)*g1
          gauss(2,l,61)=sqp(l)*g1
          gauss(3,l,61)=g6
          gauss(4,l,61)=w1
        end do
        do l=5,8
          rr=zero
          ss=zero
          if(l.eq.5) rr=-one
          if(l.eq.6) rr=one
          if(l.eq.7) ss=-one
          if(l.eq.8) ss=one
          gauss(1,l,61)=rr*g2
          gauss(2,l,61)=ss*g2
          gauss(3,l,61)=g7
          gauss(4,l,61)=w2
        end do
        do l=9,12
          gauss(1,l,61)=rqp(l-8)*g3
          gauss(2,l,61)=sqp(l-8)*g3
          gauss(3,l,61)=g4
          gauss(4,l,61)=w3
        end do
        gauss(1,13,61)=zero
        gauss(2,13,61)=zero
        gauss(3,13,61)=g5
        gauss(4,13,61)=w4
      end if
      do l=1,ngauss(61)
        rr=gauss(1,l,61)
        ss=gauss(2,l,61)
        tt=gauss(3,l,61)
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
        dq(1,1)=-fourth*((one-ss)*(-rr-ss-one)+(one-rr)*(one-ss))
        dq(2,1)= fourth*((one-ss)*( rr-ss-one)+(one+rr)*(one-ss))
        dq(3,1)= fourth*((one+ss)*( rr+ss-one)+(one+rr)*(one+ss))
        dq(4,1)=-fourth*((one+ss)*(-rr+ss-one)+(one-rr)*(one+ss))
        dq(5,1)=-rr*(one-ss)
        dq(6,1)= half*(one-ss*ss)
        dq(7,1)=-rr*(one+ss)
        dq(8,1)=-half*(one-ss*ss)
        dq(1,2)=-fourth*((one-rr)*(-rr-ss-one)+(one-rr)*(one-ss))
        dq(2,2)=-fourth*((one+rr)*( rr-ss-one)+(one+rr)*(one-ss))
        dq(3,2)= fourth*((one+rr)*( rr+ss-one)+(one+rr)*(one+ss))
        dq(4,2)= fourth*((one-rr)*(-rr+ss-one)+(one-rr)*(one+ss))
        dq(5,2)=-half*(one-rr*rr)
        dq(6,2)=-ss*(one+rr)
        dq(7,2)= half*(one-rr*rr)
        dq(8,2)=-ss*(one-rr)
        dq(9,1)=zero
        dq(9,2)=zero
        dh(1)=half*(two*tt-one)
        dh(2)=half*(two*tt+one)
        dh(3)=-two*tt
        do i=1,nen(61)
          sh(4,i,l,61)=q(iqpq(i))*h(iqph(i))
          sh(1,i,l,61)=dq(iqpq(i),1)*h(iqph(i))
          sh(2,i,l,61)=dq(iqpq(i),2)*h(iqph(i))
          sh(3,i,l,61)=q(iqpq(i))*dh(iqph(i))
        end do
      end do
      call dcopy(nshsize,sh(1,1,1,61),ione,shj(1,1,1,61),ione)
c
c...  Type 62  quadratic tetrahedron
c     r, s, and t are tetrahedral coordinates.
c
      nen(62)=10
      ngauss(62)=ione
      gauss(1,1,62)=fourth
      gauss(2,1,62)=fourth
      gauss(3,1,62)=fourth
      gauss(4,1,62)=tetvol
      if(intord.ne.2) then
        ngauss(62)=4
        g1=(five-sqrt(five))/20.0d0
        g2=(five+three*sqrt(five))/20.0d0
        do l=1,ngauss(62)
          gauss(1,l,62)=g1
          gauss(2,l,62)=g1
          gauss(3,l,62)=g1
          if(l.ne.ngauss(62)) gauss(l,l,62)=g2
        end do
      end if
      do l=1,ngauss(62)
        do i=1,nen(62)
          rr=rlt(i)*gauss(1,l,31)
          ss=slt(i)*gauss(2,l,31)
          tt=tlt(i)*gauss(3,l,31)
          uu=ulw(i)+rr+ss+tt
          sh(4,i,l,31)=uu
          sh(1,i,l,31)=rlt(i)
          sh(2,i,l,31)=slt(i)
          sh(3,i,l,31)=tlt(i)
        end do
      end do
      call dcopy(nshsize,sh(1,1,1,31),ione,shj(1,1,1,31),ione)
c
      return
      end
c
c version
c $Id: preshape-old.f,v 1.1 2004/07/06 19:19:17 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
