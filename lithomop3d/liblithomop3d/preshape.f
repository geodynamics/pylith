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
      subroutine preshape(sh,shj,gauss,intord,ietype,nen,ngauss,ierr,
     & errstrng)
c
c...program to compute information for a given element type.
c
c   The arguments are:
c     Input:
c
c       intord                  = integration order
c                                 1 = full
c                                 2 = reduced
c                                 3 = Bbar
c       ietype                  = element type
c       nen                     = number of element nodes
c       ngauss                  = number of element gauss points
c
c     Output:
c       sh(nsd+1,nen,ngauss)    = shape functions and their derivatives
c       shj(nsd+1,nen,ngauss)   = shape functions and derivatives used
c                                 for calculating Jacobian.  Unless the
c                                 element is an infinite element, this
c                                 is identical to sh.
c       gauss(nsd+1,ngauss)     = Gauss point coordinates and weights
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
      integer intord,ietype,nen,ngauss,ierr
      double precision sh(nsd+1,nen,ngauss),shj(nsd+1,nen,ngauss)
      double precision gauss(nsd+1,ngauss)
      character*(*) errstrng
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
      integer i,l,k,n,ind,nshsize,ngssize
      double precision shtmp(2160)
      double precision rr,ss,tt,drr,dss,dtt
      integer io(3)
c
cdebug      write(6,*) "Hello from preshape_f!"
c
c...  definitions
c
      ierr=izero
      ngssize=(nsd+1)*ngauss
      nshsize=ngssize*nen
c
c...  initialize arrays
c
      call fill(gauss,zero,ngssize)
      call fill(sh,zero,nshsize)
      call fill(shj,zero,nshsize)
      call fill(shtmp,zero,2160)
c
c... First type:  linear hex
c
      if(ietype.eq.1) then
        call plinhex(sh,gauss,nen,ngauss,intord)
        call dcopy(nshsize,sh,ione,shj,ione)
c
c...  Types 2-27:  linear hex + infinite boundaries
c     Use mod3 arithmetic to cover all cases
c
c******  Temporarily leave out infinite elements, since the code is not
c******  presently set up to use more than one element type.
ctemp      do n=2,27
ctemp        call dcopy(ngssize,gauss(1,1,1),ione,gauss(1,1,n),ione)
ctemp        call dcopy(nshsize,sh(1,1,1,1),ione,sh(1,1,1,n),ione)
ctemp        io(3)=(n-1)/inine
ctemp        io(2)=(n-1-io(3)*inine)/ithree
ctemp        io(1)=mod(n-1,ithree)
ctemp        do l=1,infetype(1,n)
ctemp          do i=1,infetype(2,n)
ctemp            call infellh(gauss(1,l,n),rlh(i),rr,drr,io(1))
ctemp            call infellh(gauss(2,l,n),slh(i),ss,dss,io(2))
ctemp            call infellh(gauss(3,l,n),tlh(i),tt,dtt,io(3))
ctemp            shj(4,i,l,n)=rr*ss*tt
ctemp            shj(1,i,l,n)=drr*ss*tt
ctemp            shj(2,i,l,n)=dss*rr*tt
ctemp            shj(3,i,l,n)=dtt*rr*ss
ctemp          end do
ctemp        end do
ctemp      end do
c
c...  Type 28:  linear hex with one set of collapsed nodes (wrick)
c
      else if(ietype.eq.28) then
        call plinhex(shtmp,gauss,nen+1,ngauss,intord)
        ind=1
        do l=1,ngauss
          do k=1,nen+1
            if(k.le.nen) then
              call dcopy(nsd+1,shtmp(ind),ione,sh(1,k,l),ione)
            else
              call daxpy(nsd+1,one,shtmp(ind),ione,sh(1,k-1,l),ione)
            end if
            ind=ind+nsd+1
          end do
        end do
        call dcopy(nshsize,sh,ione,shj,ione)
c
c...  Type 29:  linear hex with two sets of collapsed nodes (wedge)
c     r and s are triangular coordinates in cross-section, and t is
c     cartesian natural coordinate along the length of the wedge.
c
      else if(ietype.eq.29) then
        call plinwedge(sh,gauss,nen,ngauss,intord)
        call dcopy(nshsize,sh,ione,shj,ione)
c
c...  Type 30:  linear hex with 4 points collapsed to a point (pyramid)
c     r and s are cartesian natural coordinates on the base, and t is
c     a coordinate going from negative one at the center of the base to
c     one at the apex.
c
      else if(ietype.eq.30) then
        call plinpyr(sh,gauss,nen,ngauss,intord)
        call dcopy(nshsize,sh,ione,shj,ione)
c
c...  Type 31:  linear tetrahedron
c     r, s, and t are tetrahedral coordinates.
c     One-point integration is used in all cases.
c
      else if(ietype.eq.31) then
        call plintet(sh,gauss,nen,ngauss,intord)
        call dcopy(nshsize,sh,ione,shj,ione)
c
c... Type 32:  quadratic (20-node) hex
c
      else if(ietype.eq.32) then
        call pquadhex(sh,gauss,nen,ngauss,intord)
        call dcopy(nshsize,sh,ione,shj,ione)
c
c...  Types 33-58:  quadratic hex + infinite boundaries
c     Use mod3 arithmetic to cover all cases
c
c******  Temporarily leave out infinite elements, since the code is not
c******  presently set up to use more than one element type.
ctemp      do n=33,58
ctemp        do i=1,4
ctemp          infetype(i,n)=infetype(i,32)
ctemp        end do
ctemp        call dcopy(ngssize,gauss(1,1,32),ione,gauss(1,1,n),ione)
ctemp        call dcopy(nshsize,sh(1,1,1,32),ione,sh(1,1,1,n),ione)
ctemp        io(3)=(n-1)/inine
ctemp        io(2)=(n-1-io(3)*inine)/ithree
ctemp        io(1)=mod(n-1,ithree)
ctemp        do l=1,infetype(1,n)
ctemp          do i=1,infetype(2,n)
ctemp            call infelqh(gauss(1,l,n),rqh(i),rr,drr,io(1))
ctemp            call infelqh(gauss(2,l,n),sqh(i),ss,dss,io(2))
ctemp            call infelqh(gauss(3,l,n),tqh(i),tt,dtt,io(3))
ctemp            shj(4,i,l,n)=rr*ss*tt
ctemp            shj(1,i,l,n)=drr*ss*tt
ctemp            shj(2,i,l,n)=dss*rr*tt
ctemp            shj(3,i,l,n)=dtt*rr*ss
ctemp          end do
ctemp        end do
ctemp      end do
c
c...  Type 59:  quadratic hex with one set of collapsed nodes (wrick)
c
      else if(ietype.eq.59) then
        call pquadhex(shtmp,gauss,nen+2,ngauss,intord)
        ind=1
        do l=1,ngauss
          do k=1,nen+2
            if(k.le.7) then
              call dcopy(nsd+1,shtmp(ind),ione,sh(1,k,l),ione)
            else if(k.eq.8) then
              call daxpy(nsd+1,one,shtmp(ind),ione,sh(1,k-1,l),ione)
            else if(k.gt.13) then
              call dcopy(nsd+1,shtmp(ind),ione,sh(1,k-2,l),ione)
            else if(k.gt.8) then
              call dcopy(nsd+1,shtmp(ind),ione,sh(1,k-1,l),ione)
            end if
            ind=ind+nsd+1
          end do
        end do
        call dcopy(nshsize,sh,ione,shj,ione)
c
c...  Type 60:  quadratic hex with three sets of collapsed nodes (wedge)
c     r and s are triangular coordinates in cross-section, and t is
c     cartesian natural coordinate along the length of the wedge.
c
      else if(ietype.eq.60) then
        call pquadwedge(sh,gauss,nen,ngauss,intord)
        call dcopy(nshsize,sh,ione,shj,ione)
c
c...  Type 61:  quadratic hex with 9 points collapsed to a point
c     (pyramid).  r and s are cartesian natural coordinates on the base,
c     and t is a coordinate going from negative one at the center of the
c     base to one at the apex.
c
      else if(ietype.eq.61) then
        call pquadpyr(sh,gauss,nen,ngauss,intord)
        call dcopy(nshsize,sh,ione,shj,ione)
c
c...  Type 62  quadratic tetrahedron
c     r, s, and t are tetrahedral coordinates.
c
      else if(ietype.eq.61) then
        call pquadtet(sh,gauss,nen,ngauss,intord)
      else
        ierr=106
        errstrng="preshape"
      end if
c
      return
      end
c
c version
c $Id: preshape.f,v 1.5 2005/03/22 02:22:18 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
