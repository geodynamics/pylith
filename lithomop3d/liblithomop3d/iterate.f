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
      subroutine iterate(
     & alnz,pcg,zcg,dprev,ja,                                           ! sparse
     & bextern,btraction,bgravity,bconcforce,bprestress,bintern,bresid, ! force
     & bwork,dispvec,nforce,grav,                                       ! force
     & x,d,deld,dcur,id,iwink,wink,nsysdat,                             ! global
     & dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                     ! slip
     & nfault,dfault,tfault,                                            ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,npar,             ! elemnt
     & ielno,iside,ihistry,pres,pdir,                                   ! tractn
     & prop,mhist,infmat,infmatmod,tminmax,                             ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & histry,rtimdat,ntimdat,nvisdat,                                  ! timdat
     & rgiter,gcurr,gi,gprev,gtol,rmin,rmult,nsiter,                    ! iterate
     & skew,                                                            ! skew
     & ncodat,nunits,nprint,                                            ! ioinfo
     & getshape,bmatrix,gload_cmp,stress_cmp,stress_mat_cmp,            ! external
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
      integer ierr
      integer ja(*)
      integer id(*),iwink(*)
      integer idx(*),iwinkx(*),idslp(*),ipslp(*)
      integer nfault(*)
      integer ien(*),lm(*),lmx(*),lmf(*),infiel(*),iddmat(*)
      integer ielno(*),iside(*),ihistry(*)
      integer mhist(*),infmat(*),infmatmod(*)
      integer infetype(*)
      character errstrng*(*)
      double precision alnz(*),pcg(*),zcg(*),dprev(*)
      double precision bextern(*),btraction(*),bgravity(*),bconcforce(*)
      double precision bprestress(*),bintern(*),bresid(*),bwork(*)
      double precision dispvec(*),grav(*)
      double precision x(*),d(*),deld(*),dcur(*),wink(*)
      double precision dx(*),deldx(*),dxcur(*),winkx(*)
      double precision dfault(*),tfault(*)
      double precision s(*),stemp(*)
      double precision state(*),dstate(*),dmat(*)
      double precision pres(*),pdir(*)
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
      include "gcurr_dim.inc"
      include "gi_dim.inc"
      include "gprev_dim.inc"
      include "gtol_dim.inc"
      include "rmin_dim.inc"
      include "rmult_dim.inc"
      include "nsiter_dim.inc"
      include "ncodat_dim.inc"
      include "nunits_dim.inc"
      include "nprint_dim.inc"
c
c...  external functions
c
      external getshape,bmatrix,gload_cmp,stress_cmp,stress_mat_cmp
c
c...  intrinsic functions
c
      intrinsic mod,abs
c
c...  local variables
c
      integer i,iter
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
      include "nsiter_def.inc"
      include "ncodat_def.inc"
      include "nunits_def.inc"
      include "nprint_def.inc"
c
c...initialize convergence criteria
c
cdebug      write(6,*) "Hello from iterate_f!"
cdebug      write(6,*) "nextflag, ntractflag, ngravflag, nconcflag,"
cdebug      write(6,*) "nprestrflag, nprevdflag:"
cdebug      write(6,*) nextflag,ntractflag,ngravflag,nconcflag,
cdebug     & nprestrflag, nprevdflag
c
      fulout=.true.
      converge=.false.
      reform=.false.
      if(ireform.eq.1) reform=.true.
      skc=.false.
      do i=1,3
        gcurr(i)=10.0d0*gtol(i)
        gprev(i)=10.0d0*gtol(i)
      end do
c
c...loop over iterations
c
      do iter=1,itmaxp
        updats=((ipstrs.ne.1).or.(ipstrs.eq.1.and.nstep.eq.0.and.
     &   iter.eq.1).or.(nstep.gt.0)).and.lgdefp.ge.1
        skc=iter.gt.1.and.lgdefp.ge.1.and.(numslp.ne.0.and.
     &   (iskopt.eq.2.or.(iskopt.le.0.and.abs(iskopt).eq.nstep)))
        nittot=nittot+1
        ntimdat(6)=nittot
        reform=reform.or.(mod(iter,maxitp).eq.0)
        ireform=0
        if(reform) ireform=1
        ntimdat(9)=ireform
        used=nstep.gt.0.and.(nsol.eq.3.or.nsol.eq.4).and.iter.eq.1
        if(iter.gt.1) fulout=.false.
c
c...add pressure forces,if present, to global load vector
c
clater        if(numpr.ne.0.and.(iter.eq.1.or.updats)) call addpr(
clater     &   btot,bres,x,d,dx,tfault,histry,skew,
clater     &   ien,infin,lm,lmx,lmf,
clater     &   ielno,iside,ihistry,pres,pdir,pvec,gvec2,fulout,
clater     &   nsd,ndof,nen,nskdim,npdir,numnp,neq,nee,numrot,lastep,nhist,
clater     &   nstep,lgdefp,numel,numpr,numfn,numslp,ipstrs,idout,idebug,kto,
clater     &   kw)
c
c...add gravity body forces to global load vector
c
        if(iter.eq.1.or.updats) then
          call gload_drv(
     &     bgravity,ngravflag,grav,neq,                                 ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,                         ! elemnt
     &     prop,mhist,infmat,infmatmod,numat,npropsz,                   ! materl
     &     gauss,shj,infetype,                                          ! eltype
     &     histry,rtimdat,ntimdat,nhist,lastep,gload_cmp,               ! timdat
     &     skew,numrot,                                                 ! skew
     &     ierr,errstrng)                                               ! errcod
        end if
c
c...reform the stiffness matrix, if required
c
c**        if(reform) then
c**          if(skc) call skcomp(x,d,skew,idslp,ipslp,nsd,ndof,nskdim,
c**     &     npdim,ipstrs,numsn,numnp,nstep,lgdefp,kto)
c**          call formk(alnz,ja,
c**     &     id,idx,x,d,dx,tfault,skew,iwink,wink,iwinkx,winkx,histry,
c**     &     ien,lm,lmx,lmf,mat,infin,prop,gauss,
c**     &     dmat,stn,deps,beta,betb,iter,
c**     &     s,stemp,iddmat,
c**     &     rtimdat,ntimdat,rgiter,
c**     &     ngauss,nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,nnz,
c**     &     neq,numel,numnp,numfn,numslp,numsn,numrot,nskdim,ipstrs,
c**     &     nwink,nwinkx,nhist,lastep,idout,kto,kw)
c**        end if
c
c...if icode .eq. 1 print the stiffness matrix diagonals and stop here
c
c**        if(icode.eq.1.and.nstep.eq.0.and.iter.eq.1) then
c**          call printv(alnz,dummy,id,idx,neq,ndof,numnp,izero,idout,kw)
c**          stop
c**        else if(icode.eq.1) then
c**          stop
c**        end if
c
c...for first iteration compute residual force vector
c
c*        if(iter.eq.1) call bdiff(b,btot,bres,neq)
c
c...  compute total external load and residual force vector
c
        call bsum(bextern,btraction,bgravity,bconcforce,bprestress,
     &   bintern,bresid,nextflag,ntractflag,ngravflag,nconcflag,
     &   nprestrflag,neq)
c
c...compute the global displacement increment vector using a
c   preconditioned conjugate gradients iterative solver.  Upon
c   return the vector gvec2 contains the displacements.
c
        call pcginv(alnz,bresid,dispvec,bwork,pcg,zcg,dprev,rmin,rmult,
     &   gcurr,gprev,ja,nsiter,neq,nprevdflag,nnz,ndtot,idout,kto,kw,
     &   used)
        ntimdat(8)=ndtot
        if(nsol.eq.3.or.nsol.eq.4) then
          if(iter.eq.1) call fill(dprev,zero,neq)
          call daxpy(neq,one,dispvec,ione,dprev,ione)
        end if
c
c...for first iteration, update displacements to reflect boundary
c   conditions
c
        if(iter.eq.1) then
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
c...  zero internal force vector prior to reconputing it.
c
        call fill(bintern,zero,neq)
c
c...compute contribution to internal force vector from winkler boundary
c   conditions
c
        if(nwink.ne.0) call winklf(bintern,dispvec,iwink,wink,histry,
     &   nwink,nhist,nstep,neq,lastep)
        if(nwinkx.ne.0) call winklf(bintern,dispvec,iwinkx,winkx,histry,
     &   nwinkx,nhist,nstep,neq,lastep)
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
     &     alnz,ja,nnz,                                                 ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,iwink,wink,numnp,nwink,                                  ! global
     &     dx,iwinkx,winkx,numslp,numsn,nwinkx,                         ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     prop,mhist,infmat,infmatmod,numat,npropsz,tminmax,           ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     histry,rtimdat,rgiter,ntimdat,nhist,lastep,stress_mat_cmp,   ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else
          call stress_drv(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     prop,mhist,infmat,infmatmod,numat,npropsz,tminmax,           ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     histry,rtimdat,rgiter,ntimdat,nhist,lastep,stress_cmp,       ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        end if
c
c...compute the out-of-balance forces and convergence criteria
c
        call residu(bextern,bintern,bresid,dispvec,gtol,gi,gprev,gcurr,
     &   id,idx,neq,nextflag,numnp,iter,itmaxp,idebug,idout,kto,kw,
     &   converge)
c
c...if solution has converged, set equilibrium stresses and creep
c   strains to their current values
c
        if(converge) then
          call update_state(state,dstate,infiel,infmat,infmatmod,
     &     infetype,nstatesz,numelt,numat,ierr,errstrng)
          return
        end if
        reform=.false.
      end do
 650  format(//,"Reforming the stiffness matrix:",/)
      return
      end
c
c version
c $Id: iterate.f,v 1.7 2005/01/05 22:22:04 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
