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
      subroutine autoprestr(
     & A,                                                               ! sparse
     & bextern,btraction,bgravity,bconcforce,bintern,bresid,bwink,      ! force
     & bwinkx,dispvec,dprev,nforce,grav,                                ! force
     & x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,                      ! global
     & ibond,bond,                                                      ! bc
     & dx,deldx,dxcur,diforc,idx,iwinkx,winkx,idslp,ipslp,idhist,       ! slip
     & fault,nfault,dfault,tfault,                                      ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,           ! elemnt
     & ielno,iside,ihistry,pres,pdir,                                   ! tractn
     & prop,infmatmod,                                                  ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & histry,rtimdat,ntimdat,nvisdat,maxstp,delt,alfa,maxit,ntdinit,   ! timdat
     & lgdef,utol,ftol,etol,itmax,                                      ! timdat
     & rgiter,                                                          ! stresscmp
     & skew,                                                            ! skew
     & ncodat,nunits,nprint,istatout,nstatout,                          ! ioinfo
     & ofile,pfile,ucdroot,                                             ! files
     & autoprestrStage,iterateEvent,                                    ! PETSc logging
     & ierr,errstrng)                                                   ! errcode
c
c...  subroutine to automatically compute prestresses due to application
c     of body forces.
c
      include "implicit.inc"
c
c     PETSC includes
c
#include "finclude/petsc.h"
#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
c
c     PETSc objects
c
      PetscInt autoprestrStage
      PetscEvent iterateEvent
      Mat A
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
      integer id(*),iwink(*),iddmat(*)
      integer ibond(*)
      integer idx(*),iwinkx(*),idslp(*),ipslp(*),idhist(*)
      integer nfault(*)
      integer ien(*),lm(*),lmx(*),lmf(*),ivfamily(*)
      integer ielno(*),iside(*),ihistry(*)
      integer infmatmod(*)
      integer infetype(4)
      integer maxstp(*),maxit(*),ntdinit(*),lgdef(*),itmax(*)
      integer istatout(*),nstatout(*)
      double precision bextern(*),btraction(*),bgravity(*),bconcforce(*)
      double precision bintern(*),bresid(*),bwink(*),bwinkx(*)
      double precision dispvec(*),dprev(*),grav(*)
      double precision x(*),d(*),deld(*),dcur(*),wink(*)
      double precision bond(*)
      double precision dx(*),deldx(*),dxcur(*),diforc(*),winkx(*)
      double precision fault(*),dfault(*),tfault(*)
      double precision s(*),stemp(*)
      double precision state(*),dstate(*),state0(*),dmat(*)
      double precision  pres(*),pdir(*)
      double precision prop(*)
      double precision gauss(*),sh(*),shj(*)
      double precision histry(*),delt(*),alfa(*),utol(*),ftol(*),etol(*)
      double precision skew(*)
      character ofile*(*),pfile*(*),ucdroot*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "nforce_dim.inc"
      include "nsysdat_dim.inc"
      include "npar_dim.inc"
      include "ntimdat_dim.inc"
      include "nvisdat_dim.inc"
      include "ncodat_dim.inc"
      include "nunits_dim.inc"
      include "nprint_dim.inc"
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "gtol_dim.inc"
c
c...  external routines
c
      external bmatrixn,bmatrixb,getshapn,getshapb
      external prestr_matinit_cmp_ss,gload_cmp_ss,elas_strs_cmp_ss
      external elas_strs_mat_cmp_ss
c
c...  local variables
c
      integer igroup,naxstp,nfirst,iprestress
      integer nen,ngauss,nee,nec
      double precision time,tminmax
      logical*4 fulout,skc
c
c...  included variable definition statements
c
      include "nforce_def.inc"
      include "nsysdat_def.inc"
      include "npar_def.inc"
      include "ntimdat_def.inc"
      include "nvisdat_def.inc"
      include "ncodat_def.inc"
      include "nunits_def.inc"
      include "nprint_def.inc"
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
c
cdebug      write(6,*) "Hello from autoprestr_f!"
c
      nen=infetype(1)
      ngauss=infetype(2)
      nee=infetype(3)
      nec=infetype(4)
c
      if(idout.ne.izero) open(kw,file=ofile,status="old",
     & access="append")
      if(idsk.eq.ione) open(kp,file=pfile,status="old",access="append")
      if(idsk.eq.itwo) open(kp,file=pfile,status="old",
     & form="unformatted",access="append")
      skc=iskopt.ge.izero.and.iskopt.ne.ione.and.numslp.ne.izero
      call fill(bextern,zero,neq*nextflag)
      call fill(btraction,zero,neq*ntractflag)
      call fill(bgravity,zero,neq*ngravflag)
      call fill(bconcforce,zero,neq*nconcflag)
      call fill(bwink,zero,neq*nwinkflag)
      call fill(bwinkx,zero,neq*nwinkxflag)
      call fill(bintern,zero,neq)
      call fill(deld,zero,ndof*numnp)
      call fill(deldx,zero,ndof*numnp)
      call fill(dcur,zero,ndof*numnp)
      call fill(d,zero,ndof*numnp)
      call fill(dx,zero,ndof*numnp)
      call fill(dxcur,zero,ndof*numnp)
      call fill(state,zero,nstatesz)
      call fill(state0,zero,nstatesz0)
      call fill(dmat,zero,nddmat*ngauss*numelv)
      call fill(dstate,zero,nstatesz)
      if(numfn.ne.izero) call fill(tfault,zero,numfn*ndof)
c
      write(kto,600)
c*      call flush(kto)
      fulout=.true.
      ireform=izero
      igroup=ione
      nstep=izero
      ntimdat(1)=nstep
      naxstp=izero
      nittot=izero
      ntimdat(6)=nittot
      nrftot=izero
      ntimdat(7)=nrftot
      ndtot=izero
      ntimdat(8)=ndtot
      ntimdat(9)=ireform
      iprestress=ione
cdebug      write(6,*) "Before const:"
      call const(maxstp,delt,alfa,maxit,ntdinit,lgdef,utol,ftol,
     & etol,itmax,nintg,igroup,naxstp,nfirst,rtimdat,deltp,alfap,
     & ntimdat,nstep,maxitp,ntdinitp,lgdefp,itmaxp,gtol)
cdebug      write(6,*) "After const:"
      if(skc) then
        call skclear(idslp,skew,numsn,numnp)
        call skcomp(x,d,skew,idslp,ipslp,ipstrs,numsn,numnp,nstep,
     &   lgdefp,ierr,errstrng)
        if(ierr.ne.izero) return
      end if
c
c...  at present, no force or displacement BC are applied when computing
c     prestresses, which seems appropriate.
c
c*      call load(id,ibond,bond,dcur,deld,bconcforce,histry,deltp,numnp,
c*     & neq,nconcflag,nhist,nstep,lastep,ierr,errstrng)
cdebug      write(6,*) "After load:"
c*      if(ierr.ne.izero) return
c
c...compute current split node displacements
c
c*      if(numfn.ne.izero) then
c*        call loadf(fault,dfault,histry,deltp,nfault,nstep,numfn,nhist,
c*     &   lastep,ierr,errstrng)
c*        if(ierr.ne.izero) return
c*      end if
cdebug      write(6,*) "After loadf:"
c
c...add differential forces across internal free interfaces
c
c*      if(numdif.ne.izero) then
c*        call loadx(bconcforce,diforc,histry,idx,idhist,neq,nconcflag,
c*     &   numnp,nhist,nstep,lastep,ierr,errstrng)
c*        if(ierr.ne.izero) return
c*      end if
cdebug      write(6,*) "After loadx:"
c
c...  initialize elastic material matrices and stiffness matrix, 
c     compute forces due to applied displacements and split nodes,
c     and perform iterative solution.
c
c...  compute forces due to applied displacements and split nodes
c
      if(lgdefp.eq.izero.and.intord.ne.ithree) then
cdebug        write(6,*) "Before matinit_drv (1):"
        write(kto,650)
        call matinit_drv(
     &   A,nnz,neq,                                                     ! sparse
     &   x,d,iwink,wink,numnp,nwink,iddmat,                             ! global
     &   dx,iwinkx,winkx,numslp,numsn,nwinkx,                           ! slip
     &   tfault,numfn,                                                  ! fault
     &   s,stemp,                                                       ! stiff
     &   state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,nvfamilies,   ! elemnt
     &   numelv,nstatesz,nstatesz0,nprestrflag,ipstrs,ipauto,           ! elemnt
     &   prop,infmatmod,npropsz,tminmax,                                ! materl
     &   gauss,sh,shj,nen,ngauss,nee,                                   ! eltype
     &   histry,rtimdat,ntimdat,rgiter,nhist,lastep,                    ! timdat
     &   prestr_matinit_cmp_ss,                                         ! timdat
     &   skew,numrot,                                                   ! skew
     &   getshapn,bmatrixn,                                             ! bbar
     &   ierr,errstrng)                                                 ! errcode
cdebug        write(6,*) "After matinit_drv (1):"
c
        if(ierr.ne.izero) return
c
c*        call formdf_ss(
c*     &   bintern,neq,                                                   ! force
c*     &   x,d,dcur,numnp,iddmat,                                         ! global
c*     &   s,stemp,                                                       ! stiff
c*     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
c*     &   infmat,infmatmod,numat,                                        ! materl
c*     &   gauss,sh,shj,infetype,                                         ! eltype
c*     &   skew,numrot,                                                   ! skew
c*     &   getshapn,bmatrixn,                                             ! bbar
c*     &   ierr,errstrng)                                                 ! errcode
cdebug        write(6,*) "After formdf_ss (1):"
c
c*        if(ierr.ne.izero) return
c
c*        if(numfn.ne.izero) call formf_ss(
c*     &   bintern,neq,                                                   ! force
c*     &   x,numnp,iddmat,                                                ! global
c*     &   s,stemp,                                                       ! stiff
c*     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
c*     &   infmat,infmatmod,numat,                                        ! materl
c*     &   gauss,sh,shj,infetype,                                         ! eltype
c*     &   skew,numrot,                                                   ! skew
c*     &   nfault,dfault,tfault,numfn,                                    ! split
c*     &   getshapn,bmatrixn,                                             ! bbar
c*     &   ierr,errstrng)                                                 ! errcode
cdebug        write(6,*) "After formf_ss (1):"
c
c*        if(ierr.ne.izero) return
c
        call iterate(
     &   A,                                                             ! sparse
     &   bextern,btraction,bgravity,bconcforce,bintern,bresid,bwink,    ! force
     &   bwinkx,dispvec,dprev,nforce,grav,                              ! force
     &   x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,                    ! global
     &   dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                   ! slip
     &   nfault,dfault,tfault,                                          ! fault
     &   s,stemp,                                                       ! stiff
     &   state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,         ! elemnt
     &   ielno,iside,ihistry,pres,pdir,                                 ! tractn
     &   prop,infmatmod,tminmax,                                        ! materl
     &   gauss,sh,shj,nen,ngauss,nee,nec,                               ! eltype
     &   histry,rtimdat,ntimdat,nvisdat,iprestress,                     ! timdat
     &   rgiter,gtol,                                                   ! iterate
     &   skew,                                                          ! skew
     &   ncodat,nunits,nprint,                                          ! ioinfo
     &   getshapn,bmatrixn,gload_cmp_ss,elas_strs_cmp_ss,               ! external
     &   elas_strs_mat_cmp_ss,                                          ! external
     &   iterateEvent,                                                  ! PETSc logging
     &   ierr,errstrng)                                                 ! errcode
cdebug        write(6,*) "After iterate (1):"
c
        if(ierr.ne.izero) return
c
      else if(lgdefp.eq.izero.and.intord.eq.ithree) then
        write(kto,650)
        call matinit_drv(
     &   A,nnz,neq,                                                     ! sparse
     &   x,d,iwink,wink,numnp,nwink,iddmat,                             ! global
     &   dx,iwinkx,winkx,numslp,numsn,nwinkx,                           ! slip
     &   tfault,numfn,                                                  ! fault
     &   s,stemp,                                                       ! stiff
     &   state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,nvfamilies,   ! elemnt
     &   numelv,nstatesz,nstatesz0,nprestrflag,ipstrs,ipauto,           ! elemnt
     &   prop,infmatmod,npropsz,tminmax,                                ! materl
     &   gauss,sh,shj,nen,ngauss,nee,                                   ! eltype
     &   histry,rtimdat,ntimdat,rgiter,nhist,lastep,                    ! timdat
     &   prestr_matinit_cmp_ss,                                         ! timdat
     &   skew,numrot,                                                   ! skew
     &   getshapb,bmatrixb,                                             ! bbar
     &   ierr,errstrng)                                                 ! errcode
cdebug        write(6,*) "After matinit_drv (2):"
c
        if(ierr.ne.izero) return
c
c*        call formdf_ss(
c*     &   bintern,neq,                                                   ! force
c*     &   x,d,dcur,numnp,iddmat,                                         ! global
c*     &   s,stemp,                                                       ! stiff
c*     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
c*     &   infmat,infmatmod,numat,                                        ! materl
c*     &   gauss,sh,shj,infetype,                                         ! eltype
c*     &   skew,numrot,                                                   ! skew
c*     &   getshapb,bmatrixb,                                             ! bbar
c*     &   ierr,errstrng)                                                 ! errcode
cdebug        write(6,*) "After formdf_ss (2):"
c
c*        if(ierr.ne.izero) return
c
c*        if(numfn.ne.izero) call formf_ss(
c*     &   bintern,neq,                                                   ! force
c*     &   x,numnp,iddmat,                                                ! global
c*     &   s,stemp,                                                       ! stiff
c*     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
c*     &   infmat,infmatmod,numat,                                        ! materl
c*     &   gauss,sh,shj,infetype,                                         ! eltype
c*     &   skew,numrot,                                                   ! skew
c*     &   nfault,dfault,tfault,numfn,                                    ! split
c*     &   getshapb,bmatrixb,                                             ! bbar
c*     &   ierr,errstrng)                                                 ! errcode
cdebug        write(6,*) "After formf_ss (2):"
c
c*        if(ierr.ne.izero) return
c
        call iterate(
     &   A,                                                             ! sparse
     &   bextern,btraction,bgravity,bconcforce,bintern,bresid,bwink,    ! force
     &   bwinkx,dispvec,dprev,nforce,grav,                              ! force
     &   x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,                    ! global
     &   dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                   ! slip
     &   nfault,dfault,tfault,                                          ! fault
     &   s,stemp,                                                       ! stiff
     &   state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,         ! elemnt
     &   ielno,iside,ihistry,pres,pdir,                                 ! tractn
     &   prop,infmatmod,tminmax,                                        ! materl
     &   gauss,sh,shj,nen,ngauss,nee,nec,                               ! eltype
     &   histry,rtimdat,ntimdat,nvisdat,iprestress,                     ! timdat
     &   rgiter,gtol,                                                   ! iterate
     &   skew,                                                          ! skew
     &   ncodat,nunits,nprint,                                          ! ioinfo
     &   getshapb,bmatrixb,gload_cmp_ss,elas_strs_cmp_ss,               ! external
     &   elas_strs_mat_cmp_ss,                                          ! external
     &   iterateEvent,                                                  ! PETSc logging
     &   ierr,errstrng)                                                 ! errcode
cdebug        write(6,*) "After iterate (2):"
c
        if(ierr.ne.izero) return
c
clater      else if(lgdefp.eq.ione.and.intord.ne.ithree) then
clater        write(kto,650)
clater        call matinit_drv(
clater     &   A,nnz,neq,                                                     ! sparse
clater     &   x,d,iwink,wink,numnp,nwink,                                    ! global
clater     &   dx,iwinkx,winkx,numslp,numsn,nwinkx,                           ! slip
clater     &   tfault,numfn,                                                  ! fault
clater     &   s,stemp,                                                       ! stiff
clater     &   state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,! elemnt
clater     &   nstatesz0,ndmatsz,numelt,nconsz,nprestrflag,ipstrs,ipauto,     ! elemnt
clater     &   nstate0,                                                       ! elemnt
clater     &   prop,infmat,infmatmod,numat,npropsz,tminmax,                   ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   histry,rtimdat,ntimdat,nhist,lastep,prestr_matinit_cmp_ld,     ! timdat
clater     &   skew,numrot,                                                   ! skew
clater     &   getshapn,bmatrixn,                                             ! bbar
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
c*clater        call formdf_ld(
c*clater     &   bintern,neq,                                                   ! force
c*clater     &   x,d,dcur,numnp,iddmat,                                         ! global
c*clater     &   s,stemp,                                                       ! stiff
c*clater     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
c*clater     &   infmat,infmatmod,numat,                                        ! materl
c*clater     &   gauss,sh,shj,infetype,                                         ! eltype
c*clater     &   skew,numrot,                                                   ! skew
c*clater     &   getshapn,bmatrixn,                                             ! bbar
c*clater     &   ierr,errstrng)                                                 ! errcode
c
c*clater        if(ierr.ne.izero) return
c
c*clater        if(numfn.ne.izero) call formf_ld(
c*clater     &   bintern,neq,                                                   ! force
c*clater     &   x,numnp,iddmat,                                                ! global
c*clater     &   s,stemp,                                                       ! stiff
c*clater     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
c*clater     &   infmat,infmatmod,numat,                                        ! materl
c*clater     &   gauss,sh,shj,infetype,                                         ! eltype
c*clater     &   skew,numrot,                                                   ! skew
c*clater     &   nfault,dfault,tfault,numfn,                                    ! split
c*clater     &   getshapn,bmatrixn,                                             ! bbar
c*clater     &   ierr,errstrng)                                                 ! errcode
c
c*clater        if(ierr.ne.izero) return
c
clater        call iterate(
clater     &   A,                                                     ! sparse
clater     &   bextern,btraction,bgravity,bconcforce,bintern,bresid,bwink,    ! force
clater     &   bwinkx,dispvec,dprev,nforce,grav,                              ! force
clater     &   x,d,deld,dcur,id,iwink,wink,nsysdat,                           ! global
clater     &   dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                   ! slip
clater     &   nfault,dfault,tfault,                                          ! fault
clater     &   s,stemp,                                                       ! stiff
clater     &   state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,npar,    ! elemnt
clater     &   ielno,iside,ihistry,pres,pdir,                                 ! tractn
clater     &   prop,infmat,infmatmod,tminmax,                                 ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   histry,rtimdat,ntimdat,nvisdat,iprestress,                     ! timdat
clater     &   rgiter,gtol,                                                   ! iterate
clater     &   skew,                                                          ! skew
clater     &   ncodat,nunits,nprint,                                          ! ioinfo
clater     &   getshapn,bmatrixn,gload_cmp_ld,elas_strs_cmp_ld,               ! external
clater     &   elas_strs_mat_cmp_ld,                                          ! external
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
clater      else if(lgdefp.eq.1.and.intord.eq.ithree) then
clater        write(kto,650)
clater        call matinit_drv(
clater     &   A,nnz,neq,                                                     ! sparse
clater     &   x,d,iwink,wink,numnp,nwink,                                    ! global
clater     &   dx,iwinkx,winkx,numslp,numsn,nwinkx,                           ! slip
clater     &   tfault,numfn,                                                  ! fault
clater     &   s,stemp,                                                       ! stiff
clater     &   state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,! elemnt
clater     &   nstatesz0,ndmatsz,numelt,nconsz,nprestrflag,ipstrs,ipauto,     ! elemnt
clater     &   nstate0,                                                       ! elemnt
clater     &   prop,infmat,infmatmod,numat,npropsz,tminmax,                   ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   histry,rtimdat,ntimdat,nhist,lastep,prestr_matinit_cmp_ld,     ! timdat
clater     &   skew,numrot,                                                   ! skew
clater     &   getshapb,bmatrixb,                                             ! bbar
clater     &   ierr,errstrng)                                                 ! errcode
c
clater        if(ierr.ne.izero) return
c
c*clater        call formdf_ld(
c*clater     &   bintern,neq,                                                   ! force
c*clater     &   x,d,dcur,numnp,iddmat,                                         ! global
c*clater     &   s,stemp,                                                       ! stiff
c*clater     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
c*clater     &   infmat,infmatmod,numat,                                        ! materl
c*clater     &   gauss,sh,shj,infetype,                                         ! eltype
c*clater     &   skew,numrot,                                                   ! skew
c*clater     &   getshapb,bmatrixb,                                             ! bbar
c*clater     &   ierr,errstrng)                                                 ! errcode
c
c*clater        if(ierr.ne.izero) return
c
c*clater        if(numfn.ne.izero) call formf_ld(
c*clater     &   bintern,neq,                                                   ! force
c*clater     &   x,numnp,iddmat,                                                ! global
c*clater     &   s,stemp,                                                       ! stiff
c*clater     &   dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,           ! elemnt
c*clater     &   infmat,infmatmod,numat,                                        ! materl
c*clater     &   gauss,sh,shj,infetype,                                         ! eltype
c*clater     &   skew,numrot,                                                   ! skew
c*clater     &   nfault,dfault,tfault,numfn,                                    ! split
c*clater     &   getshapb,bmatrixb,                                             ! bbar
c*clater     &   ierr,errstrng)                                                 ! errcode
c
c*clater        if(ierr.ne.izero) return
c
clater        call iterate(
clater     &   A,                                                     ! sparse
clater     &   bextern,btraction,bgravity,bconcforce,bintern,bresid,bwink,    ! force
clater     &   bwinkx,dispvec,dprev,nforce,grav,                              ! force
clater     &   x,d,deld,dcur,id,iwink,wink,nsysdat,                           ! global
clater     &   dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,                   ! slip
clater     &   nfault,dfault,tfault,                                          ! fault
clater     &   s,stemp,                                                       ! stiff
clater     &   state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,npar,    ! elemnt
clater     &   ielno,iside,ihistry,pres,pdir,                                 ! tractn
clater     &   prop,infmat,infmatmod,tminmax,                                 ! materl
clater     &   gauss,sh,shj,infetype,                                         ! eltype
clater     &   histry,rtimdat,ntimdat,nvisdat,iprestress,                     ! timdat
clater     &   rgiter,gtol,                                                   ! iterate
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
c...print all displacements, including faulted and slippery nodes
c*    note that the following are all debugging calls, for now
c
      time=zero
      if(idsk.eq.ione) write(kp,700) nstep
      if(idsk.eq.ione) write(kp,'(e15.4)') time
      if(idsk.eq.itwo) write(kp) nstep
      if(idsk.eq.itwo) write(kp) time
      call printd(d,deld,deltp,idslp,numnp,numnp,ione,idout,idsk,kto,kw,
     & kp)
cdebug      write(6,*) "After printd (1):"
      call printf(tfault,dfault,deltp,nfault,numfn,idout,idsk,kw,kp)
cdebug      write(6,*) "After printf:"
      call printd(dx,deldx,deltp,idslp,numnp,numsn,itwo,idout,idsk,kto,
     & kw,kp)
      if(iucd.eq.ione) call write_ucd_node_vals(d,deld,tfault,dfault,
     & nfault,numfn,dx,deldx,idslp,numsn,deltp,nstep,numnp,kucd,ucdroot,
     & iprestress)
cdebug      write(6,*) "After printd (2):"
c
c...print array telling whether each slippery node is locked
c   or free for the current time step
c
      call printl(idx,iwinkx,idslp,histry,numsn,numnp,nstep,nhist,
     & nwinkx,lastep,idsk,kp)
cdebug      write(6,*) "After printl:"
c
c...print the stresses and strains
c
      call write_state_drv(
     & state,dstate,ivfamily,nvfamilies,numelv,nstatesz,                ! elemnt
     & infmatmod,                                                       ! materl
     & ngauss,                                                          ! eltype
     & delt,nstep,                                                      ! timdat
     & istatout,nstatout,                                               ! ioopts
     & idout,idsk,iucd,kw,kp,kucd,ucdroot,iprestress,                   ! ioinfo
     & ierr,errstrng)                                                   ! errcode
cdebug      write(6,*) "After write_state:"
c
      if(nintg.eq.1) then
        write(kto,800) ntimdat(6),ntimdat(7),ntimdat(8)
        if(idout.gt.0) write(kw,800) ntimdat(6),ntimdat(7),ntimdat(8)
      end if
      if(idout.ne.0) close(kw)
      close(kp)
c
c...  copy stresses into initial stress array
c
      call get_initial_stress(
     & state,state0,ivfamily,nvfamilies,numelv,nstatesz,nstatesz0,      ! elemnt
     & infmatmod,                                                       ! materl
     & nen,ngauss)                                                      ! eltype
      
c
600   format(//,'Beginning prestress computation:',/)
650   format(//,"Reforming the stiffness matrix:",/)
700   format('STEP ',i5)
800   format(/," Total number of equilibrium iterations        = ",i7,/,
     &         " Total number of stiffness matrix reformations = ",i7,/,
     &         " Total number of displacement subiterations    = ",i7)
      return
      end
c
c version
c $Id: autoprestr.F,v 1.10 2005/05/03 18:39:36 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
