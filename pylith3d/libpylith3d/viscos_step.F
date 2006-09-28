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
      subroutine viscos_step(
     & A,rhs,sol,                                                       ! sparse
     & bextern,btraction,bgravity,bconcforce,bintern,bresid,bwink,      ! force
     & bwinkx,dispvec,dprev,nforce,grav,                                ! force
     & x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,                      ! global
     & ibond,bond,                                                      ! bc
     & dx,deldx,dxcur,diforc,idx,iwinkx,winkx,idslp,ipslp,idhist,       ! slip
     & fault,nfault,dfault,tfault,                                      ! split
     & s,stemp,                                                         ! stiff
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,ielindx,   ! elemnt
     & tractionverts,tractionvals,gauss2d,sh2d,infetype2d,              ! tractn
     & prop,infmatmod,                                                  ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & histry,rtimdat,ntimdat,nvisdat,maxstp,delt,alfa,maxit,ntdinit,   ! timdat
     & lgdef,utol,ftol,etol,itmax,                                      ! timdat
     & rgiter,                                                          ! stresscmp
     & skew,                                                            ! skew
     & iprint,ncodat,nunits,nprint,istatout,nstatout,                   ! ioinfo
     & ofile,pfile,ucdroot,                                             ! files
     & viscousStage, iterateEvent,                                      ! PETSc logging
     & ntot,pyltim,indexx,jcyc,i,j,pyskc,nfirst,pynstep,time,           ! New args
     & pydeltp,lgdefp,gtol,                                             !   more new args
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
      Vec rhs, sol
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
      integer ien(*),lm(*),lmx(*),lmf(*),ivfamily(*),ielindx(*)
      integer tractionverts(*),infetype2d(*)
      integer infmatmod(*)
      integer infetype(4)
      integer maxstp(*),maxit(*),ntdinit(*),lgdef(*),itmax(*)
      integer iprint(*)
      integer istatout(*),nstatout(*)
      integer ntot, pyltim, indexx, jcyc, i, j, pyskc, nfirst, pynstep
      double precision bextern(*),btraction(*),bgravity(*),bconcforce(*)
      double precision bintern(*),bresid(*),bwink(*),bwinkx(*)
      double precision dispvec(*),dprev(*),grav(*)
      double precision x(*),d(*),deld(*),dcur(*),wink(*)
      double precision bond(*)
      double precision dx(*),deldx(*),dxcur(*),diforc(*),winkx(*)
      double precision fault(*),dfault(*),tfault(*)
      double precision s(*),stemp(*)
      double precision state(*),dstate(*),state0(*),dmat(*)
      double precision tractionvals(*),gauss2d(*),sh2d(*)
      double precision prop(*)
      double precision gauss(*),sh(*),shj(*)
      double precision histry(*),delt(*),alfa(*),utol(*),ftol(*),etol(*)
      double precision skew(*)
      double precision time, pydeltp
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
      external td_matinit_cmp_ss,gload_cmp_ss,traction_cmp_ss
      external td_strs_cmp_ss,td_strs_mat_cmp_ss
c
c...  local variables
c
cdebug      integer idb
      integer iprestress
      integer nen,ngauss,nee,nec
      double precision tminmax
      logical reform,skc,ltim,fulout,unlck,unlckf
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
c     Outside loops
c
      nen=infetype(1)
      ngauss=infetype(2)
      nee=infetype(3)
      nec=infetype(4)
      reform=.false.
      if(ireform.eq.ione) reform=.true.
      iprestress=izero
c
c     Init from inputs
c
      nstep = pynstep
      deltp = pydeltp
      skc   = pyskc.eq.1
      ltim  = pyltim.eq.1
c
c...  clear arrays at beginning of time step
c
      call fill(deld,zero,ndof*numnp)
      call fill(deldx,zero,ndof*numnp)
      call fill(bextern,zero,neq*nextflag)
      call fill(btraction,zero,neq*ntractflag)
      call fill(bconcforce,zero,neq*nconcflag)
c     call fill(bintern,zero,neq)
      if(skc) then
         call skclear(idslp,skew,numsn,numnp)
         call skcomp(x,d,skew,idslp,ipslp,ipstrs,numsn,numnp,nstep,
     &        lgdefp,ierr,errstrng)
         if(ierr.ne.izero) return
      end if
c
c...  see whether winkler forces are locked or unlocked in this step.
c
cdebug            write(6,*) "Before cklock:"
      if(nwink.ne.izero) call cklock(iwink,histry,ltim,nwink,
     &     nstep,nhist,lastep,unlck)
      if(nwinkx.ne.izero) call cklock(iwinkx,histry,ltim,nwinkx,
     &     nstep,nhist,lastep,unlckf)
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
     &     numnp,neq,nconcflag,nhist,nstep,lastep,ierr,errstrng)
cdebug            write(6,*) "After load:"
      if(ierr.ne.izero) return
c
c...  compute current split node displacements
c
      if(numfn.ne.izero) then
         call loadf(fault,dfault,histry,deltp,nfault,nstep,numfn,
     &        nhist,lastep,ierr,errstrng)
         if(ierr.ne.izero) return
      end if
c
c...  add loads from changes in differential forces across internal
c        interfaces
c
      if(numdif.ne.izero) then
         call loadx(bconcforce,diforc,histry,idx,idhist,neq,
     &        nconcflag,numnp,nhist,nstep,lastep,ierr,errstrng)
         if(ierr.ne.izero) return
      end if
c
c...  compute change in load vector if winkler forces are removed
c
      if(nwink.ne.izero.and.unlck) call unlock(bwink,bintern,
     &     iwink,histry,nstep,nwink,nhist,neq,lastep)
      if(nwinkx.ne.izero.and.unlckf) call unlock(bwinkx,bintern,
     &     iwinkx,histry,nstep,nwinkx,nhist,neq,lastep)
cdebug            write(6,*) "After unlock:"
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
     &           ipstrs,ipauto,ielindx,                                 ! elemnt
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
     &        bintern,neq,                                              ! force
     &        x,d,deld,numnp,iddmat,                                    ! global
     &        s,stemp,                                                  ! stiff
     &        dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,               ! elemnt
     &        infmatmod,                                                ! materl
     &        gauss,sh,shj,nen,ngauss,nee,                              ! eltype
     &        skew,numrot,                                              ! skew
     &        getshapn,bmatrixn,                                        ! bbar
     &        ierr,errstrng)                                            ! errcode
c     debug              write(6,*) "After formdf_ss (1):"
c     
         if(ierr.ne.izero) return
c     
         if(numfn.ne.izero) call formf_ss(
     &        bintern,neq,                                              ! force
     &        x,numnp,iddmat,                                           ! global
     &        s,stemp,                                                  ! stiff
     &        dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,               ! elemnt
     &        infmatmod,                                                ! materl
     &        gauss,sh,shj,nen,ngauss,nee,                              ! eltype
     &        skew,numrot,                                              ! skew
     &        nfault,dfault,tfault,numfn,                               ! split
     &        getshapn,bmatrixn,                                        ! bbar
     &        ierr,errstrng)                                            ! errcode
c     debug              write(6,*) "After formf_ss (1):"
c     
         if(ierr.ne.izero) return
c     
         call iterate(
     &        A,rhs,sol,                                                ! sparse
     &        bextern,btraction,bgravity,bconcforce,bintern,bresid,     ! force
     &        bwink,bwinkx,dispvec,dprev,nforce,grav,                   ! force
     &        x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,               ! global
     &        dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,              ! slip
     &        nfault,dfault,tfault,                                     ! fault
     &        s,stemp,                                                  ! stiff
     &        state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,    ! elemnt
     &        ielindx,                                                  ! elemnt
     &        tractionverts,tractionvals,gauss2d,sh2d,infetype2d,       ! tractn
     &        prop,infmatmod,tminmax,                                   ! materl
     &        gauss,sh,shj,nen,ngauss,nee,nec,                          ! eltype
     &        histry,rtimdat,ntimdat,nvisdat,iprestress,                ! timdat
     &        rgiter,gtol,                                              ! iterate
     &        skew,                                                     ! skew
     &        ncodat,nunits,nprint,                                     ! ioinfo
     &        getshapn,bmatrixn,gload_cmp_ss,traction_cmp_ss,           ! external
     &        td_strs_cmp_ss,td_strs_mat_cmp_ss,                        ! external
     &        iterateEvent,                                             ! PETSc logging
     &        ierr,errstrng)                                            ! errcode
c     debug              write(6,*) "After iterate (1):"
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
     &           ipstrs,ipauto,ielindx,                                 ! elemnt
     &           prop,infmatmod,npropsz,tminmax,                        ! materl
     &           gauss,sh,shj,nen,ngauss,nee,                           ! eltype
     &           histry,rtimdat,ntimdat,rgiter,nhist,lastep,            ! timdat
     &           td_matinit_cmp_ss,                                     ! timdat
     &           skew,numrot,                                           ! skew
     &           getshapb,bmatrixb,                                     ! bbar
     &           ierr,errstrng)                                         ! errcode
         end if
c     debug              write(6,*) "After matinit_drv (2):"
c     
         if(ierr.ne.izero) return
c     
         call formdf_ss(
     &        bintern,neq,                                              ! force
     &        x,d,deld,numnp,iddmat,                                    ! global
     &        s,stemp,                                                  ! stiff
     &        dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,               ! elemnt
     &        infmatmod,                                                ! materl
     &        gauss,sh,shj,nen,ngauss,nee,                              ! eltype
     &        skew,numrot,                                              ! skew
     &        getshapb,bmatrixb,                                        ! bbar
     &        ierr,errstrng)                                            ! errcode
c     debug              write(6,*) "After formdf_ss (2):"
c     
         if(ierr.ne.izero) return
c     
         if(numfn.ne.izero) call formf_ss(
     &        bintern,neq,                                              ! force
     &        x,numnp,iddmat,                                           ! global
     &        s,stemp,                                                  ! stiff
     &        dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,               ! elemnt
     &        infmatmod,                                                ! materl
     &        gauss,sh,shj,nen,ngauss,nee,                              ! eltype
     &        skew,numrot,                                              ! skew
     &        nfault,dfault,tfault,numfn,                               ! split
     &        getshapb,bmatrixb,                                        ! bbar
     &        ierr,errstrng)                                            ! errcode
c     debug              write(6,*) "After formf_ss (2):"
c     
         if(ierr.ne.izero) return
c     
         call iterate(
     &        A,rhs,sol,                                                ! sparse
     &        bextern,btraction,bgravity,bconcforce,bintern,bresid,     ! force
     &        bwink,bwinkx,dispvec,dprev,nforce,grav,                   ! force
     &        x,d,deld,dcur,id,iwink,wink,nsysdat,iddmat,               ! global
     &        dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,              ! slip
     &        nfault,dfault,tfault,                                     ! fault
     &        s,stemp,                                                  ! stiff
     &        state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,npar,    ! elemnt
     &        ielindx,                                                  ! elemnt
     &        tractionverts,tractionvals,gauss2d,sh2d,infetype2d,       ! tractn
     &        prop,infmatmod,tminmax,                                   ! materl
     &        gauss,sh,shj,nen,ngauss,nee,nec,                          ! eltype
     &        histry,rtimdat,ntimdat,nvisdat,iprestress,                ! timdat
     &        rgiter,gtol,                                              ! iterate
     &        skew,                                                     ! skew
     &        ncodat,nunits,nprint,                                     ! ioinfo
     &        getshapb,bmatrixb,gload_cmp_ss,traction_cmp_ss,           ! external
     &        td_strs_cmp_ss,td_strs_mat_cmp_ss,                        ! external
     &        iterateEvent,                                             ! PETSc logging
     &        ierr,errstrng)                                            ! errcode
c     debug              write(6,*) "After iterate (2):"
c     
         if(ierr.ne.izero) return
c     
c     later            else if(lgdefp.eq.ione.and.intord.ne.ithree) then
c     later              if(reform) then
c     later                write(kto,650)
c     later                call matinit_drv(
c     later     &           A,nnz,neq,                                             ! sparse
c     later     &           x,d,iwink,wink,numnp,nwink,                            ! global
c     later     &           dx,iwinkx,winkx,numslp,numsn,nwinkx,                   ! slip
c     later     &           tfault,numfn,                                          ! fault
c     later     &           s,stemp,                                               ! stiff
c     later     &           state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat, ! elemnt
c     later     &           nstatesz,nstatesz0,ndmatsz,numelt,nconsz,nprestrflag,  ! elemnt
c     later     &           ipstrs,ipauto,nstate0,                                 ! elemnt
c     later     &           prop,infmat,infmatmod,numat,npropsz,tminmax,           ! materl
c     later     &           gauss,sh,shj,infetype,                                 ! eltype
c     later     &           histry,rtimdat,ntimdat,nhist,lastep,td_matinit_cmp_ld, ! timdat
c     later     &           skew,numrot,                                           ! skew
c     later     &           getshapn,bmatrixn,                                     ! bbar
c     later     &           ierr,errstrng)                                         ! errcode
c     later              end if
c     
c     later              if(ierr.ne.izero) return
c     
c     later              call formdf_ld(
c     later     &         bintern,neq,                                             ! force
c     later     &         x,d,deld,numnp,iddmat,                                   ! global
c     later     &         s,stemp,                                                 ! stiff
c     later     &         dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,     ! elemnt
c     later     &         infmat,infmatmod,numat,                                  ! materl
c     later     &         gauss,sh,shj,infetype,                                   ! eltype
c     later     &         skew,numrot,                                             ! skew
c     later     &         getshapn,bmatrixn,                                       ! bbar
c     later     &         ierr,errstrng)                                           ! errcode
c     
c     later              if(ierr.ne.izero) return
c     
c     later              if(numfn.ne.izero) call formf_ld(
c     later     &         bintern,neq,                                             ! force
c     later     &         x,numnp,iddmat,                                          ! global
c     later     &         s,stemp,                                                 ! stiff
c     later     &         dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,     ! elemnt
c     later     &         infmat,infmatmod,numat,                                  ! materl
c     later     &         gauss,sh,shj,infetype,                                   ! eltype
c     later     &         skew,numrot,                                             ! skew
c     later     &         nfault,dfault,tfault,numfn,                              ! split
c     later     &         getshapn,bmatrixn,                                       ! bbar
c     later     &         ierr,errstrng)                                           ! errcode
c     
c     later              if(ierr.ne.izero) return
c     
c     later              call iterate(
c     later     &         A,                                               ! sparse
c     later     &         bextern,btraction,bgravity,bconcforce,bintern,bresid,    ! force
c     later     &         bwink,bwinkx,dispvec,nforce,grav,                        ! force
c     later     &         x,d,deld,dcur,id,iwink,wink,nsysdat,                     ! global
c     later     &         dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,             ! slip
c     later     &         nfault,dfault,tfault,                                    ! fault
c     later     &         s,stemp,                                                 ! stiff
c     later     &         state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,   ! elemnt
c     later     &         npar,                                                    ! elemnt
c     later     &         tractionverts,tractionvals,gauss2d,sh2d,infetype2d,      ! tractn
c     later     &         prop,infmat,infmatmod,tminmax,                           ! materl
c     later     &         gauss,sh,shj,infetype,                                   ! eltype
c     later     &         histry,rtimdat,ntimdat,nvisdat,iprestress,               ! timdat
c     later     &         rgiter,gtol,                                             ! iterate
c     later     &         skew,                                                    ! skew
c     later     &         ncodat,nunits,nprint,                                    ! ioinfo
c     later     &         getshapn,bmatrixn,gload_cmp_ld,traction_cmp_ld,          ! external
c     later     &         td_strs_cmp_ld,td_strs_mat_cmp_ld,                       ! external
c     later     &         ierr,errstrng)                                           ! errcode
c     
c     later              if(ierr.ne.izero) return
c     
c     later            else if(lgdefp.eq.ione.and.intord.eq.ithree) then
c     later              if(reform) call matinit_drv(
c     later     &         A,nnz,neq,                                               ! sparse
c     later     &         x,d,iwink,wink,numnp,nwink,                              ! global
c     later     &         dx,iwinkx,winkx,numslp,numsn,nwinkx,                     ! slip
c     later     &         tfault,numfn,                                            ! fault
c     later     &         s,stemp,                                                 ! stiff
c     later     &         state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,   ! elemnt
c     later     &         nstatesz,nstatesz0,ndmatsz,numelt,nconsz,nprestrflag,    ! elemnt
c     later     &         ipstrs,ipauto,nstate0,                                   !elemnt
c     later     &         prop,infmat,infmatmod,numat,npropsz,tminmax,             ! materl
c     later     &         gauss,sh,shj,infetype,                                   ! eltype
c     later     &         histry,rtimdat,ntimdat,nhist,lastep,td_matinit_cmp_ld,   ! timdat
c     later     &         skew,numrot,                                             ! skew
c     later     &         getshapb,bmatrixb,                                       ! bbar
c     later     &         ierr,errstrng)                                           ! errcode
c     
c     later              if(ierr.ne.izero) return
c     
c     later              call formdf_ld(
c     later     &         bintern,neq,                                             ! force
c     later     &         x,d,deld,numnp,iddmat,                                   ! global
c     later     &         s,stemp,                                                 ! stiff
c     later     &         dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,     ! elemnt
c     later     &         infmat,infmatmod,numat,                                  ! materl
c     later     &         gauss,sh,shj,infetype,                                   ! eltype
c     later     &         skew,numrot,                                             ! skew
c     later     &         getshapb,bmatrixb,                                       ! bbar
c     later     &         ierr,errstrng)                                           ! errcode
c     
c     later              if(ierr.ne.izero) return
c     
c     later              if(numfn.ne.izero) call formf_ld(
c     later     &         bintern,neq,                                             ! force
c     later     &         x,numnp,iddmat,                                          ! global
c     later     &         s,stemp,                                                 ! stiff
c     later     &         dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,     ! elemnt
c     later     &         infmat,infmatmod,numat,                                  ! materl
c     later     &         gauss,sh,shj,infetype,                                   ! eltype
c     later     &         skew,numrot,                                             ! skew
c     later     &         nfault,dfault,tfault,numfn,                              ! split
c     later     &         getshapb,bmatrixb,                                       ! bbar
c     later     &         ierr,errstrng)                                           ! errcode
c     
c     later              if(ierr.ne.izero) return
c     
c     later              call iterate(
c     later     &         A,                                               ! sparse
c     later     &         bextern,btraction,bgravity,bconcforce,bintern,bresid,    ! force
c     later     &         bwink,bwinkx,dispvec,nforce,grav,                        ! force
c     later     &         x,d,deld,dcur,id,iwink,wink,nsysdat,                     ! global
c     later     &         dx,deldx,dxcur,idx,iwinkx,winkx,idslp,ipslp,             ! slip
c     later     &         nfault,dfault,tfault,                                    ! fault
c     later     &         s,stemp,                                                 ! stiff
c     later     &         state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,   ! elemnt
c     later     &         npar,                                                    ! elemnt
c     later     &         tractionverts,tractionvals,gauss2d,sh2d,infetype2d,      ! tractn
c     later     &         prop,infmat,infmatmod,tminmax,                           ! materl
c     later     &         gauss,sh,shj,infetype,                                   ! eltype
c     later     &         histry,rtimdat,ntimdat,nvisdat,iprestress,               ! timdat
c     later     &         rgiter,gtol,                                             ! iterate
c     later     &         skew,                                                    ! skew
c     later     &         ncodat,nunits,nprint,                                    ! ioinfo
c     later     &         getshapb,bmatrixb,gload_cmp_ld,traction_cmp_ld,          ! external
c     later     &         td_strs_cmp_ld,td_strs_mat_cmp_ld,                       ! external
c     later     &         iterateEvent,                                            ! PETSc logging
c     later     &         ierr,errstrng)                                           ! errcode
c     
c     later              if(ierr.ne.izero) return
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
     &        idout,idsk,kto,kw,kp)
         call write_ucd_node_vals(d,deld,tfault,dfault,nfault,
     &        numfn,dx,deldx,idslp,numsn,deltp,ntot,numnp,kucd,iucd,
     &        ucdroot,iprestress)
         call printf(tfault,dfault,deltp,nfault,numfn,idout,
     &        idsk,kw,kp)
         call printd(dx,deldx,deltp,idslp,numnp,numsn,itwo,
     &        idout,idsk,kto,kw,kp)
         call printl(idx,iwinkx,idslp,histry,numsn,numnp,
     &        nstep,nhist,nwinkx,lastep,idsk,kp)
      end if
c
c...  print stresses and strains in all elements when requested
c
      if(fulout) then
         call write_state_drv(
     &        state,dstate,ivfamily,nvfamilies,numelv,nstatesz, ! elemnt
     &        infmatmod,        ! materl
     &        ngauss,           ! eltype
     &        deltp,ntot,       ! timdat
     &        istatout,nstatout, ! ioopts
     &        idout,idsk,iucd,kw,kp,kucd,ucdroot,iprestress, ! ioinfo
     &        ierr,errstrng)    ! errcode
      end if
      if(fulout) indexx=indexx+1
c
 650  format(//,"Reforming the stiffness matrix:",/)
 700  format('STEP ',i7)
 1000 format(////' output which follows is at time= ',1pe15.8/
     & ' step # ',i7,' in cycle # ',i7//)
 4000 format(1pe15.8)
 5000 format(///,'   Working on time ',1pe15.8,', timestep #',i7,
     & '  out of',i7)
      return
      end
c
c version
c $Id: viscos.F,v 1.15 2005/08/05 20:17:18 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
