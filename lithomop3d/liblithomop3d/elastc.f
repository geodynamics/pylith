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
      subroutine elastc(
     & alnz,pcg,zcg,ja,                                                 ! sparse
     & b,btot,bres,pvec,gvec1,gvec2,grav,                               ! force
     & x,d,deld,dprev,dcur,id,iwink,wink,nsysdat,                       ! global
     & ibond,bond,                                                      ! bc
     & dx,deldx,dxcur,diforc,idx,iwinkx,winkx,idslp,ipslp,idhist,       ! slip
     & fault,nfault,dfault,tfault,                                      ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,npar,             ! elemnt
     & ielno,iside,ihistry,pres,pdir,                                   ! tractn
     & prop,mhist,infmat,infmatmod,ismatmod,                            ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & histry,rtimdat,ntimdat,nvisdat,maxstp,delt,alfa,maxit,ntdinit,   ! timdat
     & lgdef,ibbar,utol,ftol,etol,itmax,                                ! timdat
     & rgiter,gcurr,gi,gprev,gtol,rmin,rmult,nsiter,                    ! iterate
     & skew,                                                            ! skew
     & ncodat,nunits,nprint,istatout,                                   ! ioinfo
     & ofile,pfile,                                                     ! files
     & ierr,errstrng)                                                   ! errcode
c
c...subroutine to construct and solve the elastic problem
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
      integer ja(*),id(*),iwink(*),ibond(*),idx(*),iwinkx(*)
      integer idslp(*),ipslp(*),idhist(*),nfault(*),ien(*),lm(*),lmx(*)
      integer lmf(*),infiel(*),iddmat(*),ielno(*),iside(*)
      integer ihistry(*),mhist(*),infmat(*),infmatmod(*),ismatmod(*)
      integer infetype(*),maxstp(*),maxit(*),ntdinit(*),lgdef(*)
      integer ibbar(*),itmax(*),istatout(*)
      double precision alnz(*),pcg(*),zcg(*)
      double precision b(*),btot(*),bres(*),pvec(*),gvec1(*),gvec2(*)
      double precision grav(*)
      double precision x(*),d(*),deld(*),dprev(*),dcur(*),wink(*)
      double precision bond(*)
      double precision dx(*),deldx(*),dxcur(*),diforc(*),winkx(*)
      double precision fault(*),dfault(*),tfault(*)
      double precision s(*),stemp(*)
      double precision state(*),dstate(*),dmat(*)
      double precision  pres(*),pdir(*)
      double precision prop(*)
      double precision gauss(*),sh(*),shj(*)
      double precision histry(*),delt(*),alfa(*),utol(*),ftol(*),etol(*)
      double precision skew(*)
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "nsysdat_dim.inc"
      include "npar_dim.inc"
      include "ntimdat_dim.inc"
      include "nvisdat_dim.inc"
      include "nsiter_dim.inc"
      include "ncodat_dim.inc"
      include "nunits_dim.inc"
      include "nprint_dim.inc"
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "gcurr_dim.inc"
      include "gi_dim.inc"
      include "gprev_dim.inc"
      include "gtol_dim.inc"
      include "rmin_dim.inc"
      include "rmult_dim.inc"
c
c...  external routines
c
      external bmatrixn,bmatrixb,getshapn,getshapb
      external elas_matinit_cmp_ss,gload_cmp_ss,elas_strs_cmp_ss
      external elas_strs_mat_cmp_ss
c
c...  local variables
c
      integer ii,igroup,naxstp,nfirst
      double precision time,tminmax
      logical*4 fulout,skc
c
c...  included variable definition statements
c
      include "nsysdat_def.inc"
      include "npar_def.inc"
      include "ntimdat_def.inc"
      include "nvisdat_def.inc"
      include "nsiter_def.inc"
      include "ncodat_def.inc"
      include "nunits_def.inc"
      include "nprint_def.inc"
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
c
cdebug      write(6,*) "Hello from elastc_f!"
c
      if(idout.ne.0) open(kw,file=ofile,status="old",access="append")
      if(idsk.eq.0) open(kp,file=pfile,status="old",access="append")
      if(idsk.eq.1) open(kp,file=pfile,status="old",form="unformatted",
     & access="append")
      skc=iskopt.ge.0.and.iskopt.ne.1.and.numslp.ne.0
      call fill(bres,zero,neq)
      call fill(btot,zero,neq)
      call fill(gvec1,zero,neq)
      call fill(pvec,zero,neq)
      call fill(deld,zero,ndof*numnp)
      call fill(deldx,zero,ndof*numnp)
      call fill(dcur,zero,ndof*numnp)
      call fill(d,zero,ndof*numnp)
      call fill(dx,zero,ndof*numnp)
      call fill(dxcur,zero,ndof*numnp)
      call fill(state,zero,nstr*nstatesz)
c***  Note:  This should be modified for prestresses.
      call fill(dstate,zero,nstr*nstatesz)
      if(numfn.ne.0) call fill(tfault,zero,numfn*ndof)
c
      write(kto,600)
c*      call flush(kto)
      fulout=.true.
      ireform=1
      igroup=1
      nstep=0
      ntimdat(1)=nstep
      naxstp=0
      nittot=0
      ntimdat(7)=nittot
      nrftot=0
      ntimdat(8)=nrftot
      ndtot=0
      ntimdat(9)=ndtot
      ntimdat(10)=ireform
      call const(maxstp,delt,alfa,maxit,ntdinit,lgdef,ibbar,utol,ftol,
     & etol,itmax,nintg,igroup,naxstp,nfirst,rtimdat,deltp,alfap,
     & ntimdat,nstep,maxitp,ntdinitp,lgdefp,ibbarp,itmaxp,gtol)
      if(skc) call skclear(idslp,skew,numsn,numnp)
      if(skc) call skcomp(x,d,skew,idslp,ipslp,ipstrs,numsn,numnp,nstep,
     & lgdefp,ierr,errstrng)
      if(ierr.ne.izero) return
c
c...transfer boundary conditions into global load vector btot(neq)
c   and displacement increment vector deld(ndof,numnp)
c
      call load(id,ibond,bond,dcur,deld,btot,histry,deltp,numnp,neq,
     & nhist,nstep,lastep,ierr,errstrng)
      if(ierr.ne.izero) return
c
c...compute current split node displacements
c
      if(numfn.ne.0) then
        call loadf(fault,dfault,histry,deltp,nfault,nstep,numfn,nhist,
     &   lastep,ierr,errstrng)
        if(ierr.ne.izero) return
      end if
c
c...add differential forces across internal free interfaces
c
      if(numdif.ne.0) call loadx(btot,diforc,histry,idx,idhist,neq,
     & numnp,nhist,nstep,lastep,ierr,errstrng)
      if(ierr.ne.izero) return
      ii=1
c
c...  initialize elastic material matrices and stiffness matrix, 
c     compute forces due to applied displacements and split nodes,
c     and perform iterative solution.
c
c
c************  still need to figure out how to do prestresses -- maybe
c************  with a separate program section.  Right now, it appears
c************  the best way is to compute the equivalent nodal forces and
c************  store these in a vector that gets subtracted each time.
c************  For now, ignore them and assume that they will be taken
c************  care of properly.
c*      if(nprestr.ne.0.and.ipstrs.eq.0) then
c*        call stresn(x,b,d,dx,tfault,stn,deps,beta,betb,scur,st0,dbeta,
c*     &   dbetb,skew,ien,lm,lmx,lmf,dmat,mat,prop,histry,infin,gauss,
c*     &   rtimdat,stol,iddmat,nen,numel,ndof,nsd,numnp,neq,nee,
c*     &   nstr,ngauss,nppts,ngem,nskdim,nhist,nprop,numat,numfn,numslp,
c*     &   numrot,lastep,nstep,lgdefp,ibbarp,ivisc,iplas,imhist,ipstrs,
c*     &   nprestr,nddmat,ndmat,idebug,idout,kto,kw,fulout)
c*      end if
c*************  Note that all routines need to be modified.  This shouldn't
c*************  be too hard in most cases.
c
c...  compute forces due to applied displacements and split nodes
c
c******  see about whether lgdef and ibbar should be scalars rather than vectors
c******  that vary for each time step.  If this is true, the routine names
c******  for b-bar and small strain/large deformation could be passed in at
c******  a higher level, which would simplify the code.  The only difficulty
c******  might be passing a routine name in from python/c++.
c******  Another thing to consider is that I will need to pass in some extra
c******  info that isn't currently needed for the small strain case.
c
      if(lgdefp.eq.0.and.ibbarp.eq.0) then
        call matinit_drv(
     &   alnz,ja,nnz,neq,                                               ! sparse
     &   x,d,iwink,wink,numnp,nwink,                                    ! global
     &   dx,iwinkx,winkx,numslp,numsn,nwinkx,                           ! slip
     &   tfault,numfn,                                                  ! fault
     &   s,stemp,                                                       ! stiff
     &   state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,       ! elemnt
     &   ndmatsz,numelt,nconsz,                                         ! elemnt
     &   prop,mhist,infmat,infmatmod,numat,npropsz,tminmax,             ! materl
     &   gauss,sh,shj,infetype,                                         ! eltype
     &   histry,rtimdat,ntimdat,nhist,lastep,elas_matinit_cmp_ss,       ! timdat
     &   skew,numrot,                                                   ! skew
     &   getshapn,bmatrixn,                                             ! bbar
     &   ierr,errstrng)                                                 ! errcode
c
        if(ierr.ne.izero) return
c
        call formdf_ss(
     &   b,neq,                                                         ! force
     &   x,d,dcur,numnp,                                                ! global
     &   s,stemp,                                                       ! stiff
     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
     &   infmat,infmatmod,numat,                                        ! materl
     &   gauss,sh,shj,infetype,                                         ! eltype
     &   skew,numrot,                                                   ! skew
     &   getshapn,bmatrixn,                                             ! bbar
     &   ierr,errstrng)                                                 ! errcode
c
        if(ierr.ne.izero) return
c
        if(numfn.ne.izero) call formf_ss(
     &   b,neq,                                                         ! force
     &   x,numnp,                                                       ! global
     &   s,stemp,                                                       ! stiff
     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
     &   infmat,infmatmod,numat,                                        ! materl
     &   gauss,sh,shj,infetype,                                         ! eltype
     &   skew,numrot,                                                   ! skew
     &   nfault,dfault,tfault,numfn,                                    ! split
     &   getshapn,bmatrixn,                                             ! bbar
     &   ierr,errstrng)                                                 ! errcode
c
        if(ierr.ne.izero) return
c
        call iterate(
     &   alnz,pcg,zcg,ja,                                               ! sparse
     &   b,btot,bres,pvec,gvec1,gvec2,grav,                             ! force
     &   x,d,deld,dprev,dcur,id,iwink,wink,nsysdat,                     ! global
     &   dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                   ! slip
     &   nfault,dfault,tfault,                                          ! fault
     &   s,stemp,                                                       ! stiff
     &   state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,npar,           ! elemnt
     &   ielno,iside,ihistry,pres,pdir,                                 ! tractn
     &   prop,mhist,infmat,infmatmod,                                   ! materl
     &   gauss,sh,shj,infetype,                                         ! eltype
     &   histry,rtimdat,ntimdat,nvisdat,                                ! timdat
     &   rgiter,gcurr,gi,gprev,gtol,rmin,rmult,nsiter,                  ! iterate
     &   skew,                                                          ! skew
     &   ncodat,nunits,nprint,                                          ! ioinfo
     &   getshapn,bmatrixn,gload_cmp_ss,elas_strs_cmp_ss,               ! external
     &   elas_strs_mat_cmp_ss,                                          ! external
     &   ierr,errstrng)                                                 ! errcode
c
        if(ierr.ne.izero) return
c
      else if(lgdefp.eq.0.and.ibbarp.eq.1) then
        call matinit_drv(
     &   alnz,ja,nnz,neq,                                               ! sparse
     &   x,d,iwink,wink,numnp,nwink,                                    ! global
     &   dx,iwinkx,winkx,numslp,numsn,nwinkx,                           ! slip
     &   tfault,numfn,                                                  ! fault
     &   s,stemp,                                                       ! stiff
     &   state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,       ! elemnt
     &   ndmatsz,numelt,nconsz,                                         ! elemnt
     &   prop,mhist,infmat,infmatmod,numat,npropsz,tminmax,             ! materl
     &   gauss,sh,shj,infetype,                                         ! eltype
     &   histry,rtimdat,ntimdat,nhist,lastep,elas_matinit_cmp_ss,       ! timdat
     &   skew,numrot,                                                   ! skew
     &   getshapb,bmatrixb,                                             ! bbar
     &   ierr,errstrng)                                                 ! errcode
c
        if(ierr.ne.izero) return
c
        call formdf_ss(
     &   b,neq,                                                         ! force
     &   x,d,dcur,numnp,                                                ! global
     &   s,stemp,                                                       ! stiff
     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
     &   infmat,infmatmod,numat,                                        ! materl
     &   gauss,sh,shj,infetype,                                         ! eltype
     &   skew,numrot,                                                   ! skew
     &   getshapb,bmatrixb,                                             ! bbar
     &   ierr,errstrng)                                                 ! errcode
c
        if(ierr.ne.izero) return
c
        if(numfn.ne.izero) call formf_ss(
     &   b,neq,                                                         ! force
     &   x,numnp,                                                       ! global
     &   s,stemp,                                                       ! stiff
     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
     &   infmat,infmatmod,numat,                                        ! materl
     &   gauss,sh,shj,infetype,                                         ! eltype
     &   skew,numrot,                                                   ! skew
     &   nfault,dfault,tfault,numfn,                                    ! split
     &   getshapb,bmatrixb,                                             ! bbar
     &   ierr,errstrng)                                                 ! errcode
c
        if(ierr.ne.izero) return
c
        call iterate(
     &   alnz,pcg,zcg,ja,                                               ! sparse
     &   b,btot,bres,pvec,gvec1,gvec2,grav,                             ! force
     &   x,d,deld,dprev,dcur,id,iwink,wink,nsysdat,                     ! global
     &   dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                   ! slip
     &   nfault,dfault,tfault,                                          ! fault
     &   s,stemp,                                                       ! stiff
     &   state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,npar,           ! elemnt
     &   ielno,iside,ihistry,pres,pdir,                                 ! tractn
     &   prop,mhist,infmat,infmatmod,                                   ! materl
     &   gauss,sh,shj,infetype,                                         ! eltype
     &   histry,rtimdat,ntimdat,nvisdat,                                ! timdat
     &   rgiter,gcurr,gi,gprev,gtol,rmin,rmult,nsiter,                  ! iterate
     &   skew,                                                          ! skew
     &   ncodat,nunits,nprint,                                          ! ioinfo
     &   getshapb,bmatrixb,gload_cmp_ss,elas_strs_cmp_ss,               ! external
     &   elas_strs_mat_cmp_ss,                                          ! external
     &   ierr,errstrng)                                                 ! errcode
c
        if(ierr.ne.izero) return
c
clater      else if(lgdefp.eq.1.and.ibbarp.eq.0) then
clater        call matinit_drv(
clater     &   alnz,ja,nnz,neq,                                               ! sparse
clater     &   x,d,iwink,wink,numnp,nwink,                                    ! global
clater     &   dx,iwinkx,winkx,numslp,numsn,nwinkx,                           ! slip
clater     &   tfault,numfn,                                                  ! fault
clater     &   s,stemp,                                                       ! stiff
clater     &   state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,       ! elemnt
clater     &   ndmatsz,numelt,nconsz,                                         ! elemnt
clater     &   prop,mhist,infmat,infmatmod,numat,npropsz,tminmax,             ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   histry,rtimdat,ntimdat,nhist,lastep,elas_matinit_cmp_ld,       ! timdat
clater     &   skew,numrot,                                                   ! skew
clater     &   getshapn,bmatrixn,                                             ! bbar
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
clater        call formdf_ld(
clater     &   b,neq,                                                         ! force
clater     &   x,d,dcur,numnp,                                                ! global
clater     &   s,stemp,                                                       ! stiff
clater     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
clater     &   infmat,infmatmod,numat,                                        ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   skew,numrot,                                                   ! skew
clater     &   getshapn,bmatrixn,                                             ! bbar
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
clater        if(numfn.ne.izero) call formf_ld(
clater     &   b,neq,                                                         ! force
clater     &   x,numnp,                                                       ! global
clater     &   s,stemp,                                                       ! stiff
clater     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
clater     &   infmat,infmatmod,numat,                                        ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   skew,numrot,                                                   ! skew
clater     &   nfault,dfault,tfault,numfn,                                    ! split
clater     &   getshapn,bmatrixn,                                             ! bbar
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
clater        call iterate(
clater     &   alnz,pcg,zcg,ja,                                               ! sparse
clater     &   b,btot,bres,pvec,gvec1,gvec2,grav,                             ! force
clater     &   x,d,deld,dprev,dcur,id,iwink,wink,nsysdat,                     ! global
clater     &   dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                   ! slip
clater     &   nfault,dfault,tfault,                                          ! fault
clater     &   s,stemp,                                                       ! stiff
clater     &   state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,npar,           ! elemnt
clater     &   ielno,iside,ihistry,pres,pdir,                                 ! tractn
clater     &   prop,mhist,infmat,infmatmod,                                   ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   histry,rtimdat,ntimdat,nvisdat,                                ! timdat
clater     &   rgiter,gcurr,gi,gprev,gtol,rmin,rmult,nsiter,                  ! iterate
clater     &   skew,                                                          ! skew
clater     &   ncodat,nunits,nprint,                                          ! ioinfo
clater     &   getshapn,bmatrixn,gload_cmp_ld,elas_strs_cmp_ld,               ! external
clater     &   elas_strs_mat_cmp_ld,                                          ! external
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
clater      else if(lgdefp.eq.1.and.ibbarp.eq.1) then
clater        call matinit_drv(
clater     &   alnz,ja,nnz,neq,                                               ! sparse
clater     &   x,d,iwink,wink,numnp,nwink,                                    ! global
clater     &   dx,iwinkx,winkx,numslp,numsn,nwinkx,                           ! slip
clater     &   tfault,numfn,                                                  ! fault
clater     &   s,stemp,                                                       ! stiff
clater     &   state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,       ! elemnt
clater     &   ndmatsz,numelt,nconsz,                                         ! elemnt
clater     &   prop,mhist,infmat,infmatmod,numat,npropsz,tminmax,             ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   histry,rtimdat,ntimdat,nhist,lastep,elas_matinit_cmp_ld,       ! timdat
clater     &   skew,numrot,                                                   ! skew
clater     &   getshapb,bmatrixb,                                             ! bbar
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
clater        call formdf_ld(
clater     &   b,neq,                                                         ! force
clater     &   x,d,dcur,numnp,                                                ! global
clater     &   s,stemp,                                                       ! stiff
clater     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
clater     &   infmat,infmatmod,numat,                                        ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   skew,numrot,                                                   ! skew
clater     &   getshapb,bmatrixb,                                             ! bbar
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
clater        if(numfn.ne.izero) call formf_ld(
clater     &   b,neq,                                                         ! force
clater     &   x,numnp,                                                       ! global
clater     &   s,stemp,                                                       ! stiff
clater     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
clater     &   infmat,infmatmod,numat,                                        ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   skew,numrot,                                                   ! skew
clater     &   nfault,dfault,tfault,numfn,                                    ! split
clater     &   getshapb,bmatrixb,                                             ! bbar
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
clater        call iterate(
clater     &   alnz,pcg,zcg,ja,                                               ! sparse
clater     &   b,btot,bres,pvec,gvec1,gvec2,grav,                             ! force
clater     &   x,d,deld,dprev,dcur,id,iwink,wink,nsysdat,                     ! global
clater     &   dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                   ! slip
clater     &   nfault,dfault,tfault,                                          ! fault
clater     &   s,stemp,                                                       ! stiff
clater     &   state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,npar,           ! elemnt
clater     &   ielno,iside,ihistry,pres,pdir,                                 ! tractn
clater     &   prop,mhist,infmat,infmatmod,                                   ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   histry,rtimdat,ntimdat,nvisdat,                                ! timdat
clater     &   rgiter,gcurr,gi,gprev,gtol,rmin,rmult,nsiter,                  ! iterate
clater     &   skew,                                                          ! skew
clater     &   ncodat,nunits,nprint,                                          ! ioinfo
clater     &   getshapb,bmatrixb,gload_cmp_ld,elas_strs_cmp_ld,               ! external
clater     &   elas_strs_mat_cmp_ld,                                          ! external
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
      end if
c
c...if ipstrs=1, equate prestresses with current stresses, and set
c   displacements to zero
c
c*      if(nprestr.ne.0.and.ipstrs.eq.1) then
c*        call fill(d,zero,ndof*numnp)
c*        call fill(dx,zero,ndof*numnp)
c*        call fill(deld,zero,ndof*numnp)
c*        call fill(deldx,zero,ndof*numnp)
c*        call fill(dcur,zero,ndof*numnp)
c*        call fill(dxcur,zero,ndof*numnp)
c*        call fill(deps,zero,nstr*ngauss*numel)
c*        call eqstrsql(st0,stn,x,gauss,ien,infin,nstr,ngauss,nppts,numel,
c*     &   nen,nsd,numnp,idout,kto,kw)
c*      end if
c*********Note:  in the near future, all of the calls below should be
c*********replaced by a single modular output section.
c
c...print all displacements, including faulted and slippery nodes
c
      time=zero
      if(idsk.eq.0) write(kp,700) nstep
      if(idsk.eq.0) write(kp,'(e15.4)') time
      if(idsk.eq.1) write(kp) nstep
      if(idsk.eq.1) write(kp) time
      call printd(d,deld,deltp,idslp,numnp,numnp,ione,idout,idsk,kto,kw,
     & kp)
      call printf(tfault,dfault,deltp,nfault,numfn,idout,idsk,kw,kp)
      call printd(dx,deldx,deltp,idslp,numnp,numsn,itwo,idout,idsk,kto,
     & kw,kp)
c
c...print array telling whether each slippery node is locked
c   or free for the current time step
c
      call printl(idx,iwinkx,idslp,histry,numsn,numnp,nstep,nhist,
     & nwinkx,lastep,idsk,kp)
c
c...print the stresses and strains
c
      call write_state(
     & state,dstate,infiel,nstatesz,numelt,                             ! elemnt
     & infmat,infmatmod,ismatmod,numat,                                 ! materl
     & infetype,                                                        ! eltype
     & delt,nstep,                                                      ! timdat
     & istatout,                                                        ! ioopts
     & idout,idsk,kw,kp)                                                ! ioinfo
c
      if(nintg.eq.1) then
        write(kto,800) ntimdat(7),ntimdat(8),ntimdat(9)
        if(idout.gt.0) write(kw,800) ntimdat(7),ntimdat(8),ntimdat(9)
      end if
      if(idout.ne.0) close(kw)
      close(kp)
c
600   format(//,'Beginning elastic solution:',/)
700   format('STEP ',i5)
800   format(/," Total number of equilibrium iterations        = ",i7,/,
     &         " Total number of stiffness matrix reformations = ",i7,/,
     &         " Total number of displacement subiterations    = ",i7)
      return
      end
c
c version
c $Id: elastc.f,v 1.3 2004/07/13 16:11:16 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
