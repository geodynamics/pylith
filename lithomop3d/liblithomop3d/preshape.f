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
      subroutine preshape(sh,shj,gauss,infetype,intord)
c
c...program to compute the shape functions in natural coordinates.
c   The program also assigns integration points and weights.
c   All shape functions are stored in a single array.
c
c   The arguments are:
c     Input:
c
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
c       infetype(4,netypes)                = element type info array
c         infetype(1,i) = number of gauss points for element type
c         infetype(2,i) = number of nodes for element type
c         infetype(3,i) = number of element coordinate directions for
c                         element type
c         infetype(4,i) = number of element degrees of freedom for
c                         element type
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
      integer intord
      integer infetype(4,netypes)
      double precision sh(nsd+1,nenmax,ngaussmax,netypes)
      double precision shj(nsd+1,nenmax,ngaussmax,netypes)
      double precision gauss(nsd+1,ngaussmax,netypes)
c
c...  local constants
c
      double precision rlh(8),slh(8),tlh(8)
      data rlh/-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0/
      data slh/-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0/
      data tlh/ 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0/
c
      double precision rqh(20),sqh(20),tqh(20)
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
c...  intrinsic functions
c
      intrinsic mod
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i,l,n,nshsize,ngssize
      double precision rr,ss,tt,drr,dss,dtt
      integer io(3)
c
cdebug      write(6,*) "Hello from preshape_f!"
c
c...  definitions
c
      nshsize=(nsd+1)*nenmax*ngaussmax
      ngssize=(nsd+1)*ngaussmax
c
c...  initialize arrays
c
      call ifill(infetype,izero,ifour*netypes)
      call fill(gauss,zero,ngssize*netypes)
      call fill(sh,zero,nshsize*netypes)
      call fill(shj,zero,nshsize*netypes)
c
c... First type:  linear hex
c
      call plinhex(sh(1,1,1,1),shj(1,1,1,1),gauss(1,1,1),infetype(1,1),
     & intord)
c
c...  Types 2-27:  linear hex + infinite boundaries
c     Use mod3 arithmetic to cover all cases
c
      do n=2,27
        do i=1,4
          infetype(i,n)=infetype(i,1)
        end do
        call dcopy(ngssize,gauss(1,1,1),ione,gauss(1,1,n),ione)
        call dcopy(nshsize,sh(1,1,1,1),ione,sh(1,1,1,n),ione)
        io(3)=(n-1)/inine
        io(2)=(n-1-io(3)*inine)/ithree
        io(1)=mod(n-1,ithree)
        do l=1,infetype(1,n)
          do i=1,infetype(2,n)
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
      infetype(1,28)=infetype(1,1)
      infetype(2,28)=iseven
      infetype(3,28)=nsd*infetype(2,28)
      infetype(4,28)=ndof*infetype(2,28)
      call dcopy(ngssize,gauss(1,1,1),ione,gauss(1,1,28),ione)
      call dcopy(nshsize,sh(1,1,1,1),ione,sh(1,1,1,28),ione)
      do l=1,infetype(1,28)
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
      call plinwedge(sh(1,1,1,29),shj(1,1,1,29),gauss(1,1,29),
     & infetype(1,29),intord)
c
c...  Type 30:  linear hex with 4 points collapsed to a point (pyramid)
c     r and s are cartesian natural coordinates on the base, and t is
c     a coordinate going from negative one at the center of the base to
c     one at the apex.
c
      call plinpyr(sh(1,1,1,30),shj(1,1,1,30),gauss(1,1,30),
     & infetype(1,30),intord)
c
c...  Type 31:  linear tetrahedron
c     r, s, and t are tetrahedral coordinates.
c     One-point integration is used in all cases.
c
      call plintet(sh(1,1,1,31),shj(1,1,1,31),gauss(1,1,31),
     & infetype(1,31),intord)
c
c... Type 32:  quadratic (20-node) hex
c
      call pquadhex(sh(1,1,1,32),shj(1,1,1,32),gauss(1,1,32),
     & infetype(1,32),intord)
c
c...  Types 33-58:  quadratic hex + infinite boundaries
c     Use mod3 arithmetic to cover all cases
c
      do n=33,58
        do i=1,4
          infetype(i,n)=infetype(i,32)
        end do
        call dcopy(ngssize,gauss(1,1,32),ione,gauss(1,1,n),ione)
        call dcopy(nshsize,sh(1,1,1,32),ione,sh(1,1,1,n),ione)
        io(3)=(n-1)/inine
        io(2)=(n-1-io(3)*inine)/ithree
        io(1)=mod(n-1,ithree)
        do l=1,infetype(1,n)
          do i=1,infetype(2,n)
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
      infetype(1,59)=infetype(1,32)
      infetype(2,59)=18
      infetype(3,59)=nsd*infetype(2,59)
      infetype(4,59)=ndof*infetype(2,59)
      call dcopy(ngssize,gauss(1,1,32),ione,gauss(1,1,59),ione)
c*      call dcopy(nshsize,sh(1,1,1,32),ione,sh(1,1,1,59),ione)
c*      do l=1,infetype(1,59)
c*        do i=1,4
c*          sh(i,7,l,59)=sh(i,7,l,59)+sh(i,8,l,59)+sh(i,15,l,59)
c*          sh(i,8,l,59)=zero
c*          sh(i,15,l,59)=zero
c*        end do
c*      end do
      do n=1,infetype(2,59)
        do l=1,infetype(1,59)
          do i=1,4
            sh(i,n,l,59)=sh(i,n,l,32)
            if(n.eq.7) sh(i,n,l,59)=sh(i,7,l,32)+sh(i,8,l,32)+
     &       sh(i,15,l,32)
            if(n.gt.7) sh(i,n,l,59)=sh(i,n+1,l,32)
            if(n.gt.13) sh(i,n,l,59)=sh(i,n+2,l,32)
          end do
        end do
      end do
      call dcopy(nshsize,sh(1,1,1,59),ione,shj(1,1,1,59),ione)
c
c...  Type 60:  quadratic hex with three sets of collapsed nodes (wedge)
c     r and s are triangular coordinates in cross-section, and t is
c     cartesian natural coordinate along the length of the wedge.
c
      call pquadwedge(sh(1,1,1,60),shj(1,1,1,60),gauss(1,1,60),
     & infetype(1,60),intord)
c
c...  Type 61:  quadratic hex with 9 points collapsed to a point
c     (pyramid).  r and s are cartesian natural coordinates on the base,
c     and t is a coordinate going from negative one at the center of the
c     base to one at the apex.
c
      call pquadpyr(sh(1,1,1,61),shj(1,1,1,61),gauss(1,1,61),
     & infetype(1,61),intord)
c
c...  Type 62  quadratic tetrahedron
c     r, s, and t are tetrahedral coordinates.
c
      call pquadtet(sh(1,1,1,62),shj(1,1,1,62),gauss(1,1,62),
     & infetype(1,62),intord)
c
      return
      end
c
c version
c $Id: preshape.f,v 1.4 2004/08/13 16:42:48 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
