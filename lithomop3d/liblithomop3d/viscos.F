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
      subroutine viscos(
     & A,                                                               ! sparse
     & bextern,btraction,bgravity,bconcforce,bintern,bresid,dispvec,    ! force
     & dprev,nforce,grav,                                               ! force
     & x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,                      ! global
     & ibond,bond,                                                      ! bc
     & dx,deldx,dxcur,diforc,idx,iwinkx,winkx,idslp,ipslp,idhist,       ! slip
     & fault,nfault,dfault,tfault,                                      ! split
     & s,stemp,                                                         ! stiff
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,           ! elemnt
     & ielno,iside,ihistry,pres,pdir,                                   ! tractn
     & prop,infmatmod,                                                  ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & histry,rtimdat,ntimdat,nvisdat,maxstp,delt,alfa,maxit,ntdinit,   ! timdat
     & lgdef,utol,ftol,etol,itmax,                                      ! timdat
     & rgiter,                                                          ! stresscmp
     & skew,                                                            ! skew
     & iprint,ncodat,nunits,nprint,istatout,nstatout,                   ! ioinfo
     & ofile,pfile,ucdroot,                                             ! files
     & viscousStage, iterateEvent,                                      ! PETSc logging
     & ierr,errstrng)                                                   ! errcode
c
c...subroutine to solve the time dependent problem and perform the
c   time stepping
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
      PetscInt viscousStage
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
      integer iprint(*)
      integer istatout(*),nstatout(*)
      double precision bextern(*),btraction(*),bgravity(*),bconcforce(*)
      double precision bintern(*),bresid(*),dispvec(*),dprev(*),grav(*)
      double precision x(*),d(*),deld(*),dcur(*),wink(*)
      double precision bond(*)
      double precision dx(*),deldx(*),dxcur(*),diforc(*),winkx(*)
      double precision fault(*),dfault(*),tfault(*)
      double precision s(*),stemp(*)
      double precision state(*),dstate(*),state0(*),dmat(*)
      double precision pres(*),pdir(*)
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
c...  intrinsic functions
c
      intrinsic abs,mod
c
c...  external routines
c
      external bmatrixn,bmatrixb,getshapn,getshapb
      external td_matinit_cmp_ss,gload_cmp_ss,td_strs_cmp_ss
      external td_strs_mat_cmp_ss
c
c...  local variables
c
cdebug      integer idb
      integer indexx,ntot,jcyc,nfirst,naxstp,i,j,iprestress
      integer nen,ngauss,nee,nec
      double precision time,tminmax
      logical ltim,fulout,unlck,unlckf,skc,reform
c
c...  included variable definitions
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
c...  open output files for appending, if necessary
c
cdebug      write(6,*) "Hello from viscos_f!"
c
      call PetscLogStagePush(viscousStage, ierr)
c
      nen=infetype(1)
      ngauss=infetype(2)
      nee=infetype(3)
      nec=infetype(4)
c
      if(idout.gt.ione) open(kw,file=ofile,status="old",access="append")
      if(idsk.eq.ione) open(kp,file=pfile,status="old",access="append")
      if(idsk.eq.itwo) open(kp,file=pfile,status="old",
     & form="unformatted",access="append")
c
c...  signal user that viscous computation is begun
c
      write(kto,2000)
c
c...  loop over complete cycles
c
      reform=.false.
      if(ireform.eq.ione) reform=.true.
      indexx=ione
      ntot=izero
      iprestress=izero
      do jcyc=1,ncycle
        if(ncycle.gt.ione) write(kto,2001) jcyc
        nfirst=izero
        naxstp=izero
        nstep=izero
        time=zero
c
c...  loop over time step groups
c
        do i=2,nintg
c
c...  define constants to control stepping in current group
c
cdebug          write(6,*) "Before const:"
          call const(maxstp,delt,alfa,maxit,ntdinit,lgdef,utol,
     &     ftol,etol,itmax,nintg,i,naxstp,nfirst,rtimdat,deltp,alfap,
     &     ntimdat,nstep,maxitp,ntdinitp,lgdefp,itmaxp,gtol)
cdebug          write(6,*) "After const:"
cdebug          write(6,*) nintg,i,naxstp,nfirst,deltp,alfap,maxitp,maxitcp,
cdebug     &     lgdefp,itmaxp,(gtol(idb),idb=1,3)
          ltim=.true.
c
c...  loop over time steps in current group
c
          do j=nfirst,naxstp
            ntot=ntot+ione
            nstep=nstep+ione
            ntimdat(1)=nstep
            time=time+deltp
            skc=(numslp.ne.0.and.(iskopt.eq.2.or.(iskopt.le.0.and.
     &       abs(iskopt).eq.nstep)))
c
c...  clear arrays at beginning of time step
c
            call fill(deld,zero,ndof*numnp)
            call fill(deldx,zero,ndof*numnp)
            call fill(bextern,zero,neq*nextflag)
            call fill(btraction,zero,neq*ntractflag)
            call fill(bconcforce,zero,neq*nconcflag)
c*            call fill(bintern,zero,neq)
            if(skc) then
              call skclear(idslp,skew,numsn,numnp)
              call skcomp(x,d,skew,idslp,ipslp,ipstrs,numsn,numnp,nstep,
     &         lgdefp,ierr,errstrng)
              if(ierr.ne.izero) return
            end if
c
c...  see whether winkler forces are locked or unlocked in this step.
c
cdebug            write(6,*) "Before cklock:"
            if(nwink.ne.izero) call cklock(iwink,histry,ltim,nwink,
     &       nstep,nhist,lastep,unlck)
            if(nwinkx.ne.izero) call cklock(iwinkx,histry,ltim,nwinkx,
     &       nstep,nhist,lastep,unlckf)
cdebug            write(6,*) "After cklock:"
c
c...  test for reform and refactor interval, whether full output
c     occurs in this step
c
            ireform=izero
            if(ntdinitp.eq.izero) then
              reform=.false.
            else if(ntdinitp.lt.izero) then
              reform=.false.
              if(j.eq.nfirst) reform=.true.
            else
              reform=(mod(j,ntdinitp).eq.izero)
            end if
            reform=reform.or.ltim
            if(reform) ireform=ione
            ntimdat(9)=ireform
cdebug            write(6,*) "ntdinitp,reform,ireform,ltim:",
cdebug     &       ntdinitp,reform,ireform,ltim
            fulout=.false.
            if(ntot.eq.iprint(indexx)) fulout=.true.
c
            if(idout.gt.ione) write(kw,1000) time,ntot,jcyc
C***********************************
            if(fulout.and.idsk.eq.ione) write(kp,700) ntot
C***********************************
            if(fulout.and.idsk.eq.ione) write(kp,4000) time
            if(fulout.and.idsk.eq.itwo) write(kp) ntot
            if(fulout.and.idsk.eq.itwo) write(kp) time
            write(kto,5000) time,ntot,lastep*ncycle
            call flush(kto)
c*            call flush(kw)
c*            call flush(kp)
c
c...  apply boundary conditions
c
cdebug            write(6,*) "Before load:"
            call load(id,ibond,bond,d,deld,bconcforce,histry,deltp,
     &       numnp,neq,nconcflag,nhist,nstep,lastep,ierr,errstrng)
cdebug            write(6,*) "After load:"
            if(ierr.ne.izero) return
c
c...  compute current split node displacements
c
            if(numfn.ne.izero) then
              call loadf(fault,dfault,histry,deltp,nfault,nstep,numfn,
     &         nhist,lastep,ierr,errstrng)
              if(ierr.ne.izero) return
            end if
c
c...  add loads from changes in differential forces across internal
c        interfaces
c
            if(numdif.ne.izero) then
              call loadx(bconcforce,diforc,histry,idx,idhist,neq,
     &         nconcflag,numnp,nhist,nstep,lastep,ierr,errstrng)
              if(ierr.ne.izero) return
            end if
c
c...  compute change in load vector if winkler forces are removed
c
c****  I believe that this is no longer necessary, since Winkler forces
c****  are added to the internal force vector at every iteration.
c
c*            call fill(zcg,zero,neq)
c*            if(nwink.ne.izero.and.unlck) call unlock(zcg,btot,id,iwink,
c*     &       idhist,ibond,bond,histry,nstep,numnp,nwink,nhist,neq,
c*     &       numdif,lastep,ione)
c*            if(nwinkx.ne.izero.and.unlckf) call unlock(zcg,btot,idx,
c*     &       iwinkx,idhist,ibond,diforc,histry,nstep,numnp,nwinkx,nhist,
c*     &       neq,numdif,lastep,itwo)
cdebug            write(6,*) "After unlock:"
c*	    call daxpy(neq,one,zcg,ione,btot,ione)
c
c...  reform time-dependent material and stiffness matrices if
c     requested, compute forces due to applied displacements and split
c     nodes, and perform iterative solution.
c
            if(lgdefp.eq.izero.and.intord.ne.ithree) then
cdebug              write(6,*) "Before matinit_drv (1):"
cdebug              write(6,*) "reform:",reform
              if(reform) then
                write(kto,650)
                call matinit_drv(
     &           A,nnz,neq,                                             ! sparse
     &           x,d,iwink,wink,numnp,nwink,iddmat,                     ! global
     &           dx,iwinkx,winkx,numslp,numsn,nwinkx,                   ! slip
     &           tfault,numfn,                                          ! fault
     &           s,stemp,                                               ! stiff
     &           state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,      ! elemnt
     &           nvfamilies,numelv,nstatesz,nstatesz0,nprestrflag,      ! elemnt
     &           ipstrs,ipauto,                                         ! elemnt
     &           prop,infmatmod,npropsz,tminmax,                        ! materl
     &           gauss,sh,shj,nen,ngauss,nee,                           ! eltype
     &           histry,rtimdat,ntimdat,rgiter,nhist,lastep,            ! timdat
     &           td_matinit_cmp_ss,                                     ! timdat
     &           skew,numrot,                                           ! skew
     &           getshapn,bmatrixn,                                     ! bbar
     &           ierr,errstrng)                                         ! errcode
              end if
cdebug              write(6,*) "After matinit_drv (1):"
c
              if(ierr.ne.izero) return
c
              call formdf_ss(
     &         bintern,neq,                                             ! force
     &         x,d,dcur,numnp,iddmat,                                   ! global
     &         s,stemp,                                                 ! stiff
     &         dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,              ! elemnt
     &         infmatmod,                                               ! materl
     &         gauss,sh,shj,nen,ngauss,nee,                             ! eltype
     &         skew,numrot,                                             ! skew
     &         getshapn,bmatrixn,                                       ! bbar
     &         ierr,errstrng)                                           ! errcode
cdebug              write(6,*) "After formdf_ss (1):"
c
              if(ierr.ne.izero) return
c
              if(numfn.ne.izero) call formf_ss(
     &         bintern,neq,                                             ! force
     &         x,numnp,iddmat,                                          ! global
     &         s,stemp,                                                 ! stiff
     &         dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,              ! elemnt
     &         infmatmod,                                               ! materl
     &         gauss,sh,shj,nen,ngauss,nee,                             ! eltype
     &         skew,numrot,                                             ! skew
     &         nfault,dfault,tfault,numfn,                              ! split
     &         getshapn,bmatrixn,                                       ! bbar
     &         ierr,errstrng)                                           ! errcode
cdebug              write(6,*) "After formf_ss (1):"
c
              if(ierr.ne.izero) return
c
              call iterate(
     &         A,                                                       ! sparse
     &         bextern,btraction,bgravity,bconcforce,bintern,bresid,    ! force
     &         dispvec,dprev,nforce,grav,                               ! force
     &         x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,              ! global
     &         dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,             ! slip
     &         nfault,dfault,tfault,                                    ! fault
     &         s,stemp,                                                 ! stiff
     &         state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,   ! elemnt
     &         ielno,iside,ihistry,pres,pdir,                           ! tractn
     &         prop,infmatmod,tminmax,                                  ! materl
     &         gauss,sh,shj,nen,ngauss,nee,nec,                         ! eltype
     &         histry,rtimdat,ntimdat,nvisdat,iprestress,               ! timdat
     &         rgiter,gtol,                                             ! iterate
     &         skew,                                                    ! skew
     &         ncodat,nunits,nprint,                                    ! ioinfo
     &         getshapn,bmatrixn,gload_cmp_ss,td_strs_cmp_ss,           ! external
     &         td_strs_mat_cmp_ss,                                      ! external
     &         iterateEvent,                                            ! PETSc logging
     &         ierr,errstrng)                                           ! errcode
cdebug              write(6,*) "After iterate (1):"
c
              if(ierr.ne.izero) return
c
            else if(lgdefp.eq.izero.and.intord.eq.ithree) then
cdebug              write(6,*) "Before matinit_drv (2):"
cdebug              write(6,*) "reform:",reform
              if(reform) then
                write(kto,650)
                call matinit_drv(
     &           A,nnz,neq,                                             ! sparse
     &           x,d,iwink,wink,numnp,nwink,iddmat,                     ! global
     &           dx,iwinkx,winkx,numslp,numsn,nwinkx,                   ! slip
     &           tfault,numfn,                                          ! fault
     &           s,stemp,                                               ! stiff
     &           state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,      ! elemnt
     &           nvfamilies,numelv,nstatesz,nstatesz0,nprestrflag,      ! elemnt
     &           ipstrs,ipauto,                                         ! elemnt
     &           prop,infmatmod,npropsz,tminmax,                        ! materl
     &           gauss,sh,shj,nen,ngauss,nee,                           ! eltype
     &           histry,rtimdat,ntimdat,rgiter,nhist,lastep,            ! timdat
     &           td_matinit_cmp_ss,                                     ! timdat
     &           skew,numrot,                                           ! skew
     &           getshapb,bmatrixb,                                     ! bbar
     &           ierr,errstrng)                                         ! errcode
              end if
cdebug              write(6,*) "After matinit_drv (2):"
c
              if(ierr.ne.izero) return
c
              call formdf_ss(
     &         bintern,neq,                                             ! force
     &         x,d,dcur,numnp,iddmat,                                   ! global
     &         s,stemp,                                                 ! stiff
     &         dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,              ! elemnt
     &         infmatmod,                                               ! materl
     &         gauss,sh,shj,nen,ngauss,nee,                             ! eltype
     &         skew,numrot,                                             ! skew
     &         getshapb,bmatrixb,                                       ! bbar
     &         ierr,errstrng)                                           ! errcode
cdebug              write(6,*) "After formdf_ss (2):"
c
              if(ierr.ne.izero) return
c
              if(numfn.ne.izero) call formf_ss(
     &         bintern,neq,                                             ! force
     &         x,numnp,iddmat,                                          ! global
     &         s,stemp,                                                 ! stiff
     &         dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,              ! elemnt
     &         infmatmod,                                               ! materl
     &         gauss,sh,shj,nen,ngauss,nee,                             ! eltype
     &         skew,numrot,                                             ! skew
     &         nfault,dfault,tfault,numfn,                              ! split
     &         getshapb,bmatrixb,                                       ! bbar
     &         ierr,errstrng)                                           ! errcode
cdebug              write(6,*) "After formf_ss (2):"
c
              if(ierr.ne.izero) return
c
              call iterate(
     &         A,                                                       ! sparse
     &         bextern,btraction,bgravity,bconcforce,bintern,bresid,    ! force
     &         dispvec,dprev,nforce,grav,                               ! force
     &         x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,              ! global
     &         dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,             ! slip
     &         nfault,dfault,tfault,                                    ! fault
     &         s,stemp,                                                 ! stiff
     &         state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,   ! elemnt
     &         ielno,iside,ihistry,pres,pdir,                           ! tractn
     &         prop,infmatmod,tminmax,                                  ! materl
     &         gauss,sh,shj,nen,ngauss,nee,nec,                         ! eltype
     &         histry,rtimdat,ntimdat,nvisdat,iprestress,               ! timdat
     &         rgiter,gtol,                                             ! iterate
     &         skew,                                                    ! skew
     &         ncodat,nunits,nprint,                                    ! ioinfo
     &         getshapb,bmatrixb,gload_cmp_ss,td_strs_cmp_ss,           ! external
     &         td_strs_mat_cmp_ss,                                      ! external
     &         iterateEvent,                                            ! PETSc logging
     &         ierr,errstrng)                                           ! errcode
cdebug              write(6,*) "After iterate (2):"
c
              if(ierr.ne.izero) return
c
clater            else if(lgdefp.eq.ione.and.intord.ne.ithree) then
clater              if(reform) then
clater                write(kto,650)
clater                call matinit_drv(
clater     &           A,nnz,neq,                                             ! sparse
clater     &           x,d,iwink,wink,numnp,nwink,                            ! global
clater     &           dx,iwinkx,winkx,numslp,numsn,nwinkx,                   ! slip
clater     &           tfault,numfn,                                          ! fault
clater     &           s,stemp,                                               ! stiff
clater     &           state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat, ! elemnt
clater     &           nstatesz,nstatesz0,ndmatsz,numelt,nconsz,nprestrflag,  ! elemnt
clater     &           ipstrs,ipauto,nstate0,                                 ! elemnt
clater     &           prop,infmat,infmatmod,numat,npropsz,tminmax,           ! materl
clater     &           gauss,sh,shj,infetype,                                 ! eltype
clater     &           histry,rtimdat,ntimdat,nhist,lastep,td_matinit_cmp_ld, ! timdat
clater     &           skew,numrot,                                           ! skew
clater     &           getshapn,bmatrixn,                                     ! bbar
clater     &           ierr,errstrng)                                         ! errcode
clater              end if
c
clater              if(ierr.ne.izero) return
c
clater              call formdf_ld(
clater     &         bintern,neq,                                             ! force
clater     &         x,d,dcur,numnp,iddmat,                                   ! global
clater     &         s,stemp,                                                 ! stiff
clater     &         dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,     ! elemnt
clater     &         infmat,infmatmod,numat,                                  ! materl
clater     &         gauss,sh,shj,infetype,                                   ! eltype
clater     &         skew,numrot,                                             ! skew
clater     &         getshapn,bmatrixn,                                       ! bbar
clater     &         ierr,errstrng)                                           ! errcode
c
clater              if(ierr.ne.izero) return
c
clater              if(numfn.ne.izero) call formf_ld(
clater     &         bintern,neq,                                             ! force
clater     &         x,numnp,iddmat,                                          ! global
clater     &         s,stemp,                                                 ! stiff
clater     &         dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,     ! elemnt
clater     &         infmat,infmatmod,numat,                                  ! materl
clater     &         gauss,sh,shj,infetype,                                   ! eltype
clater     &         skew,numrot,                                             ! skew
clater     &         nfault,dfault,tfault,numfn,                              ! split
clater     &         getshapn,bmatrixn,                                       ! bbar
clater     &         ierr,errstrng)                                           ! errcode
c
clater              if(ierr.ne.izero) return
c
clater              call iterate(
clater     &         A,                                               ! sparse
clater     &         bextern,btraction,bgravity,bconcforce,bintern,bresid,    ! force
clater     &         dispvec,nforce,grav,                                     ! force
clater     &         x,d,deld,dcur,id,iwink,wink,nsysdat,                     ! global
clater     &         dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,             ! slip
clater     &         nfault,dfault,tfault,                                    ! fault
clater     &         s,stemp,                                                 ! stiff
clater     &         state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,   ! elemnt
clater     &         npar,                                                    ! elemnt
clater     &         ielno,iside,ihistry,pres,pdir,                           ! tractn
clater     &         prop,infmat,infmatmod,tminmax,                           ! materl
clater     &         gauss,sh,shj,infetype,                                   ! eltype
clater     &         histry,rtimdat,ntimdat,nvisdat,iprestress,               ! timdat
clater     &         rgiter,gtol,                                             ! iterate
clater     &         skew,                                                    ! skew
clater     &         ncodat,nunits,nprint,                                    ! ioinfo
clater     &         getshapn,bmatrixn,gload_cmp_ld,td_strs_cmp_ld,           ! external
clater     &         td_strs_mat_cmp_ld,                                      ! external
clater     &         ierr,errstrng)                                           ! errcode
c
clater              if(ierr.ne.izero) return
c
clater            else if(lgdefp.eq.ione.and.intord.eq.ithree) then
clater              if(reform) call matinit_drv(
clater     &         A,nnz,neq,                                               ! sparse
clater     &         x,d,iwink,wink,numnp,nwink,                              ! global
clater     &         dx,iwinkx,winkx,numslp,numsn,nwinkx,                     ! slip
clater     &         tfault,numfn,                                            ! fault
clater     &         s,stemp,                                                 ! stiff
clater     &         state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,   ! elemnt
clater     &         nstatesz,nstatesz0,ndmatsz,numelt,nconsz,nprestrflag,    ! elemnt
clater     &         ipstrs,ipauto,nstate0,                                   !elemnt
clater     &         prop,infmat,infmatmod,numat,npropsz,tminmax,             ! materl
clater     &         gauss,sh,shj,infetype,                                   ! eltype
clater     &         histry,rtimdat,ntimdat,nhist,lastep,td_matinit_cmp_ld,   ! timdat
clater     &         skew,numrot,                                             ! skew
clater     &         getshapb,bmatrixb,                                       ! bbar
clater     &         ierr,errstrng)                                           ! errcode
c
clater              if(ierr.ne.izero) return
c
clater              call formdf_ld(
clater     &         bintern,neq,                                             ! force
clater     &         x,d,dcur,numnp,iddmat,                                   ! global
clater     &         s,stemp,                                                 ! stiff
clater     &         dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,     ! elemnt
clater     &         infmat,infmatmod,numat,                                  ! materl
clater     &         gauss,sh,shj,infetype,                                   ! eltype
clater     &         skew,numrot,                                             ! skew
clater     &         getshapb,bmatrixb,                                       ! bbar
clater     &         ierr,errstrng)                                           ! errcode
c
clater              if(ierr.ne.izero) return
c
clater              if(numfn.ne.izero) call formf_ld(
clater     &         bintern,neq,                                             ! force
clater     &         x,numnp,iddmat,                                          ! global
clater     &         s,stemp,                                                 ! stiff
clater     &         dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,     ! elemnt
clater     &         infmat,infmatmod,numat,                                  ! materl
clater     &         gauss,sh,shj,infetype,                                   ! eltype
clater     &         skew,numrot,                                             ! skew
clater     &         nfault,dfault,tfault,numfn,                              ! split
clater     &         getshapb,bmatrixb,                                       ! bbar
clater     &         ierr,errstrng)                                           ! errcode
c
clater              if(ierr.ne.izero) return
c
clater              call iterate(
clater     &         A,                                               ! sparse
clater     &         bextern,btraction,bgravity,bconcforce,bintern,bresid,    ! force
clater     &         dispvec,nforce,grav,                                     ! force
clater     &         x,d,deld,dcur,id,iwink,wink,nsysdat,                     ! global
clater     &         dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,             ! slip
clater     &         nfault,dfault,tfault,                                    ! fault
clater     &         s,stemp,                                                 ! stiff
clater     &         state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,   ! elemnt
clater     &         npar,                                                    ! elemnt
clater     &         ielno,iside,ihistry,pres,pdir,                           ! tractn
clater     &         prop,infmat,infmatmod,tminmax,                           ! materl
clater     &         gauss,sh,shj,infetype,                                   ! eltype
clater     &         histry,rtimdat,ntimdat,nvisdat,iprestress,               ! timdat
clater     &         rgiter,gtol,                                             ! iterate
clater     &         skew,                                                    ! skew
clater     &         ncodat,nunits,nprint,                                    ! ioinfo
clater     &         getshapb,bmatrixb,gload_cmp_ld,td_strs_cmp_ld,           ! external
clater     &         td_strs_mat_cmp_ld,                                      ! external
clater     &         iterateEvent,                                            ! PETSc logging
clater     &         ierr,errstrng)                                           ! errcode
c
clater              if(ierr.ne.izero) return
c
            end if
c
c...  transfer current displacement increments into displacement vector
c     if a nonzero starting guess is desired.
c
            if(nprevdflag.eq.1) call idisp(deld,dprev,id,numnp,neq)
c
c...  print displacements at all nodes when requested.
c
            if(fulout) then
              call printd(d,deld,deltp,idslp,numnp,numnp,ione,
     &         idout,idsk,kto,kw,kp)
              call write_ucd_node_vals(d,deld,tfault,dfault,nfault,
     &         numfn,dx,deldx,idslp,numsn,deltp,nstep,numnp,kucd,
     &         ucdroot,iprestress)
              call printf(tfault,dfault,deltp,nfault,numfn,idout,
     &         idsk,kw,kp)
              call printd(dx,deldx,deltp,idslp,numnp,numsn,itwo,
     &         idout,idsk,kto,kw,kp)
              call printl(idx,iwinkx,idslp,histry,numsn,numnp,
     &         nstep,nhist,nwinkx,lastep,idsk,kp)
            end if
c
c...  print stresses and strains in all elements when requested
c
            if(fulout) then
              call write_state_drv(
     &         state,dstate,ivfamily,nvfamilies,numelv,nstatesz,        ! elemnt
     &         infmatmod,                                               ! materl
     &         ngauss,                                                  ! eltype
     &         deltp,nstep,                                             ! timdat
     &         istatout,nstatout,                                       ! ioopts
     &         idout,idsk,iucd,kw,kp,kucd,ucdroot,iprestress,           ! ioinfo
     &         ierr,errstrng)                                           ! errcode
            end if
            ltim=.false.
            if(fulout) indexx=indexx+1
            if(indexx.gt.icontr) indexx=icontr
          end do
        end do
      end do
      write(kto,800) ntimdat(6),ntimdat(7),ntimdat(8)
      if(idout.gt.ione) write(kw,800) ntimdat(6),ntimdat(7),ntimdat(8)
      if(idout.gt.ione) close(kw)
      close(kp)
c
 650  format(//,"Reforming the stiffness matrix:",/)
 700  format('STEP ',i7)
 800  format(/," Total number of equilibrium iterations        = ",i7,/,
     &         " Total number of stiffness matrix reformations = ",i7,/,
     &         " Total number of displacement subiterations    = ",i7)
 1000 format(////' output which follows is at time= ',1pe15.8/
     & ' step # ',i7,' in cycle # ',i7//)
 2000 format(/'Time-dependent solution is begun:'/)
 2001 format('     working on cycle ',i7)
 4000 format(1pe15.8)
 5000 format(///,'   Working on time ',1pe15.8,', timestep #',i7,
     & '  out of',i7)
      call PetscLogStagePop(ierr)
      return
      end
c
c version
c $Id: viscos.F,v 1.10 2005/04/16 00:38:59 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
