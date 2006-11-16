c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine iterate(
     & A,rhs,sol,                                                       ! sparse
     & bextern,btraction,bgravity,bconcforce,bintern,bresid,bwink,      ! force
     & bwinkx,dispvec,dprev,nforce,grav,                                ! force
     & x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,                      ! global
     & dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                     ! slip
     & nfault,dfault,tfault,                                            ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,ielindx,   ! elemnt
     & tractionverts,tractionvals,gauss2d,sh2d,infetype2d,              ! tractn
     & prop,infmatmod,tminmax,                                          ! materl
     & gauss,sh,shj,nen,ngauss,nee,nec,                                 ! eltype
     & histry,rtimdat,ntimdat,nvisdat,iprestress,                       ! timdat
     & rgiter,gtol,                                                     ! iterate
     & skew,                                                            ! skew
     & ncodat,nunits,nprint,                                            ! ioinfo
     & getshape,bmatrix,gload_cmp,traction_cmp,                         ! external
     & stress_cmp,stress_mat_cmp,                                       ! external
     & iterateEvent,                                                    ! PETSc logging
     & ierr,errstrng)                                                   ! errcode
c
c...subroutine to loop over iterations within each time step
c
c       This is a program segment of tecton that performs an iterative
c       solution until both the force tolerance and the energy tolerance
c       are achieved, or until the maximum number of iterations has
c       been exceeded.  In the case of large deformations, pressure BC
c       and body forces are re-evaluated over the new geometry.
c       This routine is called by both the elastic and viscoelastic
c       portions of the code.  The appropriate routines to use for
c       b-bar/no b-bar, small/large strain, and elastic/time-dependent
c       calculations are determined by the external routine names that
c       are passed in.
c
      include "implicit.inc"
c
c     PETSC includes
c
#include "finclude/petsc.h"
#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscviewer.h"
c
c     PETSc objects
c
      PetscEvent iterateEvent
      PetscTruth flag
      PetscInt useSolver
      Vec rhs, sol
      Vec locRhs, locSol
      Mat A
      KSP ksp
      PetscInt row, nz
      PetscViewer rhsCheck, solCheck
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
      include "rconsts.inc"

c
c...  subroutine arguments
c
      integer nen,ngauss,nee,nec,iprestress,ierr
      integer id(*),iwink(*),iddmat(*)
      integer idx(*),iwinkx(*),idslp(*),ipslp(*)
      integer nfault(*)
      integer ien(*),lm(*),lmx(*),lmf(*),ivfamily(*)
      integer tractionverts(*),infetype2d(*)
      integer ielindx(*)
      integer infmatmod(*)
      character errstrng*(*)
      double precision bextern(*),btraction(*),bgravity(*),bconcforce(*)
      double precision bintern(*),bresid(*),bwink(*),bwinkx(*)
      double precision dispvec(*),dprev(*),grav(*)
      double precision x(*),d(*),deld(*),dcur(*),wink(*)
      double precision dx(*),deldx(*),dxcur(*),winkx(*)
      double precision dfault(*),tfault(*)
      double precision s(*),stemp(*)
      double precision state(*),dstate(*),state0(*),dmat(*)
      double precision tractionvals(*),gauss2d(*),sh2d(*)
      double precision prop(*),tminmax
      double precision gauss(*),sh(*),shj(*)
      double precision histry(*)
      double precision skew(*)
c
c...  included dimension and type statements
c
      include "nforce_dim.inc"
      include "nsysdat_dim.inc"
      include "npar_dim.inc"
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
      include "nvisdat_dim.inc"
      include "rgiter_dim.inc"
      include "gtol_dim.inc"
      include "ncodat_dim.inc"
      include "nunits_dim.inc"
      include "nprint_dim.inc"
c
c...  external functions
c
      external getshape,bmatrix,gload_cmp,traction_cmp
      external stress_cmp,stress_mat_cmp
c
c...  intrinsic functions
c
      intrinsic mod,abs
c
c...  local variables
c
      integer i,iter,nsnodes,nsgauss
      double precision gi(3),gprev(3),gcurr(3)
      logical fulout,converge,updats,skc,used,reform
c
c...  included variable definitions
c
      include "nforce_def.inc"
      include "nsysdat_def.inc"
      include "npar_def.inc"
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
      include "nvisdat_def.inc"
      include "rgiter_def.inc"
      include "ncodat_def.inc"
      include "nunits_def.inc"
      include "nprint_def.inc"
c
c...initialize convergence criteria
c
cdebug      write(6,*) "Hello from iterate_f!"
c
      call PetscLogEventBegin(iterateEvent, ierr)
      nsnodes=infetype2d(1)
      nsgauss=infetype2d(2)
      fulout=.true.
      converge=.false.
      reform=.false.
      if(ireform.eq.1) reform=.true.
      skc=.false.
      do i=1,3
        gcurr(i)=10.0d0*gtol(i)
        gprev(i)=10.0d0*gtol(i)
      end do
cdebug      write(6,*) "reform,ireform:",reform,ireform
c
c...loop over iterations
c
      do iter=1,itmaxp
        updats=iter.eq.1.or.lgdefp.ge.1
        skc=iter.gt.1.and.lgdefp.ge.1.and.(numslp.ne.0.and.
     &   (iskopt.eq.2.or.(iskopt.le.0.and.abs(iskopt).eq.nstep)))
        nittot=nittot+1
        ntimdat(6)=nittot
        reform=reform.or.(mod(iter,maxitp).eq.0)
        ireform=0
        if(reform) ireform=1
        ntimdat(9)=ireform
        used=nstep.gt.0.and.nprevdflag.eq.1.and.iter.eq.1
        if(iter.gt.1) fulout=.false.
c
c...add traction forces,if present, to global load vector
c
        if(numtractions.ne.0.and.updats.and.iprestress.eq.0) then
          call traction_drv(
     &     btraction,ntractflag,neq,                                    ! force
     &     x,d,id,numnp,                                                ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! split
     &     gauss2d,sh2d,nsnodes,nsgauss,traction_cmp,                   ! eltype
     &     tractionverts,tractionvals,numtractions,                     ! traction
     &     skew,numrot,                                                 ! skew
     &     ierr,errstrng)                                               ! errcode
        end if
c
c...add gravity body forces to global load vector
c
        if(updats) then
          call gload_drv(
cpetsc     &     rhs,                                                         ! PETSc
     &     bgravity,ngravflag,grav,neq,                                 ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     ien,lm,lmx,lmf,ivfamily,nvfamilies,numelv,                   ! elemnt
     &     prop,infmatmod,npropsz,                                      ! materl
     &     gauss,shj,nen,ngauss,nee,                                    ! eltype
     &     histry,rtimdat,ntimdat,nhist,lastep,gload_cmp,               ! timdat
     &     skew,numrot,                                                 ! skew
     &     ierr,errstrng)                                               ! errcod
        end if
c
c...  compute total external load and residual force vector
c
        call bsum(bextern,btraction,bgravity,bconcforce,bintern,bresid,
     &   nextflag,ntractflag,ngravflag,nconcflag,neq)
c
c...compute the global displacement increment vector using a
c   preconditioned conjugate gradients iterative solver.  Upon
c   return the vector dispvec contains the displacements.
c
ccc PETSc BLOCK
        call VecCreateSeqWithArray(MPI_COMM_SELF, neq, bresid, locRhs,
     &                               ierr)
        call VecZeroEntries(rhs, ierr)
        call assembleVectorComplete(rhs, locRhs, ADD_VALUES, ierr)
        call VecDestroy(locRhs, ierr)
cdebug        call VecView(rhs, PETSC_VIEWER_STDOUT_WORLD, ierr)
        call KSPCreate(MPI_COMM_WORLD, ksp, ierr)
        call KSPSetfromOptions(ksp, ierr)
        if (used) then
          call KSPSetInitialGuessNonzero(ksp, 1, ierr)
          call dcopy(neq, dprev, 1, dispvec, 1)
        end if
        call KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)
        call KSPSolve(ksp, rhs, sol, ierr)
        call KSPDestroy(ksp, ierr)
cdebug        call VecView(sol, PETSC_VIEWER_STDOUT_WORLD, ierr)
        call VecCreateSeqWithArray(MPI_COMM_SELF, neq, dispvec, locSol,
     &   ierr)
        call restrictVector(sol, locSol, INSERT_VALUES, ierr)
        call VecDestroy(locSol, ierr)
c
c...for first iteration, update displacements to reflect boundary
c   conditions
c
        if(iter.eq.1.and.iprestress.eq.0) then
          if(numfn.ne.0.and.numrot.ne.0) call rsplit(nfault,dfault,
     &     skew,numfn,numnp)
          if(numfn.ne.0) call daxpy(ndof*numfn,one,dfault,ione,tfault,
     &     ione)
          if(nstep.eq.0) then
            if(numrot.ne.0) call rdisp(dcur,skew,numnp)
            call dcopy(ndof*numnp,dcur,ione,d,ione)
            call dcopy(ndof*numnp,dcur,ione,deld,ione)
          else
            if(numrot.ne.0) call rdisp(deld,skew,numnp)
            call daxpy(ndof*numnp,one,deld,ione,d,ione)
          end if
        end if
c
c...localize the displacement increments in dcur(ndof,numnp) and
c   dxcur(ndof,numnp)
c
        call fill(dcur,zero,ndof*numnp)
        call fill(dxcur,zero,ndof*numnp)
        call disp(dispvec,dcur,id,numnp,neq)
        if(numslp.ne.0) call disp(dispvec,dxcur,idx,numnp,neq)
c
c...rotate skewed coordinates to global system
c
        if(numrot.ne.0) then
          call rdisp(dcur,skew,numnp)
          if(numslp.ne.0) call rdisp(dxcur,skew,numnp)
        end if
c
c...  zero internal force vector prior to recomputing it.
c
        call fill(bintern,zero,neq)
c
c...compute contribution to internal force vector from winkler boundary
c   conditions
c
        if(nwink.ne.0) call winklf(bwink,bintern,dispvec,iwink,wink,
     &   histry,nwink,nhist,nstep,neq,lastep)
        if(nwinkx.ne.0) call winklf(bwinkx,bintern,dispvec,iwinkx,winkx,
     &   histry,nwinkx,nhist,nstep,neq,lastep)
c
c...update the total displacement and the displacement increment
c   for this step
c
        call daxpy(ndof*numnp,one,dcur,ione,d,ione)
        call daxpy(ndof*numnp,one,dcur,ione,deld,ione)
        if(numslp.ne.0) call daxpy(ndof*numnp,one,dxcur,ione,dx,ione)
        if(numslp.ne.0) call daxpy(ndof*numnp,one,dxcur,ione,deldx,ione)
c
c...  integrate the stresses and compute the equivalent nodal loads,
c     updating the material and stiffness matrices if requested.
c
        if(reform) then
          write(kto,650)
          call stress_mat_drv(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,iwink,wink,numnp,nwink,iddmat,                           ! global
     &     dx,iwinkx,winkx,numslp,numsn,nwinkx,                         ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,nvfamilies, ! elemnt
     &     numelv,nstatesz,nstatesz0,nprestrflag,ipstrs,ipauto,ielindx, ! elemnt
     &     prop,infmatmod,npropsz,tminmax,                              ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     histry,rtimdat,rgiter,ntimdat,nhist,lastep,stress_mat_cmp,   ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
          write(6,*) "tminmax:  ",tminmax
        else
          call stress_drv(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,nvfamilies, ! elemnt
     &     numelv,nstatesz,nstatesz0,nprestrflag,ipstrs,ipauto,         ! elemnt
     &     prop,infmatmod,npropsz,tminmax,                              ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     histry,rtimdat,rgiter,ntimdat,nhist,lastep,stress_cmp,       ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
          write(6,*) "tminmax:  ",tminmax
        end if
        if(ierr.ne.izero) return
c
c...compute the out-of-balance forces and convergence criteria
c
        call residu(rhs,bextern,bintern,bresid,dispvec,
     &   gtol,gi,gprev,gcurr,
     &   id,idx,neq,nextflag,numnp,iter,itmaxp,idebug,idout,kto,kw,
     &   converge)
c
c...if solution has converged, set equilibrium stresses and creep
c   strains to their current values
c
        if(converge) then
          call update_state_drv(state,dstate,ivfamily,infmatmod,
     &     nvfamilies,numelv,nstatesz,ngauss,ierr,errstrng)
          call PetscLogEventEnd(iterateEvent, ierr)
          return
        end if
        reform=.false.
      end do
 650  format(//,"Reforming the stiffness matrix:",/)
      call PetscLogEventEnd(iterateEvent, ierr)
      return
      end
c
c version
c $Id: iterate.F,v 1.11 2005/05/03 18:40:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
