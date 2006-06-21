c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
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
      subroutine residu(rhs,bextern,bintern,bresid,dispvec,gtol,gi,
     & gprev,
     & gcurr,id,idx,neq,nextflag,numnp,iter,itmaxp,idebug,idout,kto,kw,
     & converge)
c
c...computes the deviation from overall force and energy balance
c
      include "implicit.inc"
c     PETSC includes
c
#include "finclude/petsc.h"
#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
c
c     PETSc objects
c
      Vec rhs
      Vec locRhs, locSol
      PetscErrorCode ierr
c
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
      integer iout
      parameter(iout=1)
      double precision alpha
      parameter(alpha=-1.0d0)
c
c...  subroutine arguments
c
      integer neq,nextflag,numnp,iter,itmaxp,idebug,idout,kto,kw
      integer id(ndof,numnp),idx(ndof,numnp)
      double precision bextern(nextflag*neq),bintern(neq),bresid(neq)
      double precision dispvec(neq),gtol(3),gi(3),gprev(3),gcurr(3)
      logical converge
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
      intrinsic sqrt,max
c
c...  user-defined functions
c
      double precision dnrm2
      external dnrm2
c
c...  local variables
c
cdebug      integer idb
      integer i
      double precision ftmp,en,rnum
      double precision tmp(3),acc(3)
      logical debug,div(3)
c
      call flush(6)
cdebug      write(6,*) "Hello from residu_f!"
c
      tmp(3)=zero
      debug=(idebug.eq.1).and.(idout.gt.1)
      if(debug) call printv(bresid,bextern,id,idx,neq,numnp,iout,idout,
     & kw)
      call VecCreateSeqWithArray(MPI_COMM_SELF, neq, bresid,  locRhs,
     &  ierr)
      call VecCreateSeqWithArray(MPI_COMM_SELF, neq, dispvec, locSol,
     &  ierr)
      call VecZeroEntries(rhs, ierr)
      call assembleVectorComplete(rhs, locRhs, ADD_VALUES, ierr)
cdebug      write(6,*) "Initial bresid:",(bresid(idb),idb=1,neq)
cdebug      write(6,*) "Initial RHS:"
cdebug      call flush(6)
cdebug      call VecView(rhs, PETSC_VIEWER_STDOUT_WORLD, ierr)
      call VecNorm(rhs, NORM_2, ftmp, ierr)
cpetsc      ftmp=dnrm2(neq,bresid,ione)
      do i=1,neq
        en=dispvec(i)*bresid(i)
        tmp(3)=tmp(3)+en*en
      end do
      if(nextflag.eq.izero) then
        call fill(bresid,zero,neq)
      else
        call dcopy(neq,bextern,ione,bresid,ione)
      end if
      call daxpy(neq,alpha,bintern,ione,bresid,ione)
      call VecZeroEntries(rhs, ierr)
      call assembleVectorComplete(rhs, locSol, ADD_VALUES, ierr)
      call VecNorm(rhs, NORM_2, tmp(1), ierr)
cpetsc tmp(1)=dnrm2(neq,dispvec,ione)
      call VecZeroEntries(rhs, ierr)
      call assembleVectorComplete(rhs, locRhs, ADD_VALUES, ierr)
cdebug      write(6,*) "Final bresid:",(bresid(idb),idb=1,neq)
cdebug      write(6,*) "Final RHS:"
cdebug      call flush(6)
cdebug      call VecView(rhs, PETSC_VIEWER_STDOUT_WORLD, ierr)
      call VecNorm(rhs, NORM_2, tmp(2), ierr)
cpetsc      tmp(2)=dnrm2(neq,bresid,ione)
      tmp(3)=sqrt(tmp(3))
      call VecDestroy(locRhs, ierr)
      call VecDestroy(locSol, ierr)
      if(iter.eq.1) then
        gi(1)=tmp(1)
        gi(2)=max(ftmp,tmp(2))
        gi(3)=tmp(3)
        do i=1,3
          rnum=one
          if(gi(i).ne.zero) rnum=gi(i)
          gcurr(i)=rnum
          gprev(i)=rnum
        end do
      else if(iter.eq.2) then
        do i=1,3
          if(tmp(i).gt.gi(i)) then
            gi(i)=tmp(i)
            gcurr(i)=tmp(i)
            gprev(i)=tmp(i)
          end if
        end do
      end if
      do i=1,3
        if(gi(i).eq.zero) gi(i)=tmp(i)
        acc(i)=zero
        if(gi(i).ne.zero) acc(i)=tmp(i)/gi(i)
        div(i)=tmp(i).gt.gcurr(i).and.gcurr(i).gt.gprev(i)
        gprev(i)=gcurr(i)
        gcurr(i)=tmp(i)
      end do
      if(div(1).and.div(2).and.div(3)) then
        converge=.true.
        if(idout.gt.1) write(kw,800) iter
        write(kto,800) iter
      else if(acc(1).le.gtol(1).and.acc(2).le.gtol(2).and.acc(3).le.
     & gtol(3)) then
        converge=.true.
        if(idout.gt.1) write(kw,810) iter
        write(kto,810) iter
      else if(iter.eq.itmaxp) then
        converge=.true.
        if(idout.gt.1) write(kw,820) iter
        write(kto,820) iter
      end if
      if(converge) then
        if(idout.gt.1) write(kw,900) gcurr(1),gcurr(2),gcurr(3),gi(1),
     &   gi(2),gi(3),acc(1),acc(2),acc(3)
        write(kto,900) gcurr(1),gcurr(2),gcurr(3),gi(1),gi(2),gi(3),
     &   acc(1),acc(2),acc(3)
      end if
800   format(/,
     & '   WARNING!  Apparent divergence in iteration #',i5,'!',//)
810   format(/,
     & '   Equilibrium achieved in',i5,' iterations',//)
820   format(/,
     & '   WARNING!  Solution failed to converge in',i5,' iterations!',
     & //)
900   format(23x,'Displacement',9x,'Force',13x,'Energy',/,3x,
     & '   Final Norm      ',1pe15.8,3x,1pe15.8,3x,1pe15.8,/,3x,
     & '   Initial Norm    ',1pe15.8,3x,1pe15.8,3x,1pe15.8,/,3x,
     & '   Final/Initial   ',1pe15.8,3x,1pe15.8,3x,1pe15.8,//)
      return
      end
c
c version
c $Id: residu.f,v 1.4 2005/03/29 21:43:32 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
