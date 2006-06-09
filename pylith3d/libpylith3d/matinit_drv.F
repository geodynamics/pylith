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
      subroutine matinit_drv(
     & A,nnz,neq,                                                       ! sparse
     & x,d,iwink,wink,numnp,nwink,iddmat,                               ! global
     & dx,iwinkx,winkx,numslp,numsn,nwinkx,                             ! slip
     & tfault,numfn,                                                    ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,nvfamilies,     ! elemnt
     & numelv,nstatesz,nstatesz0,nprestrflag,ipstrs,ipauto,ielindx,     ! elemnt
     & prop,infmatmod,npropsz,tminmax,                                  ! materl
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & histry,rtimdat,ntimdat,rgiter,nhist,lastep,matinit_cmp,          ! timdat
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...  driver subroutine to form the initial material matrices and
c     stiffness matrix
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
      integer nnz,neq,numnp,nwink,numslp,numsn,nwinkx,numfn,nvfamilies
      integer numelv,nstatesz,nstatesz0,nprestrflag,ipstrs,ipauto
      integer npropsz,nen,ngauss,nee,nhist,lastep,numrot,ierr
      integer iwink(2,nwink),iddmat(nstr,nstr),iwinkx(2,nwinkx)
      integer ien(nen,numelv),ielindx(numelv)
      integer lm(ndof*nen,numelv),lmx(ndof*nen,numelv),lmf(nen,numelv)
      integer ivfamily(5,nvfamilies)
      integer infmatmod(6,nmatmodmax)
      character errstrng*(*)
      double precision x(nsd,numnp),d(ndof,numnp)
      double precision wink(nwink),dx(ndof,numnp),winkx(nwinkx)
      double precision tfault(ndof,numfn)
      double precision s(neemax*neemax),stemp(neemax*neemax)
      double precision state(nstatesz),dstate(nstatesz)
      double precision state0(nstatesz0)
      double precision dmat(nddmat*ngauss,numelv),prop(npropsz),tminmax
      double precision gauss(nsd+1,ngauss)
      double precision sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss)
      double precision histry(nhist,lastep+1),skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
      include "rgiter_dim.inc"
c
c...  external routines
c
      include "elas_matinit_ext.inc"
      include "td_matinit_ext.inc"
      include "prestr_matinit_ext.inc"
      external matinit_cmp,getshape,bmatrix
c
c...  local variables
c
      integer ielg,ifam,nelfamily,matmodel,indstate,indstate0,indprop
      integer nstate,nprop,nstate0,n0states,imat
      double precision ptmp(100)
      logical matchg
cdebug      integer idb,jdb,kdb
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
      include "rgiter_def.inc"
c
cdebug      write(6,*) "Hello from matinit_drv_f!"
c
      tminmax=big
      nrftot=ntimdat(7)
      nrftot=nrftot+1
      ntimdat(7)=nrftot
      ielg=ione
c
c...  Zero all entries in stiffness matrix
c
      call MatZeroEntries(A, ierr)
c
c...  loop over element families and then select appropriate material
c     model routine
c
      do ifam=1,nvfamilies
        nelfamily=ivfamily(1,ifam)
        matmodel=ivfamily(2,ifam)
        indstate=ivfamily(3,ifam)
        indstate0=ivfamily(4,ifam)
        indprop=ivfamily(5,ifam)
        nstate=infmatmod(2,matmodel)
        nprop=infmatmod(3,matmodel)
        nstate0=infmatmod(6,matmodel)
        n0states=ione
        if(ipstrs.ne.izero) n0states=nelfamily
        imat=ifam
        matchg=.false.
        call dcopy(nprop,prop(indprop),ione,ptmp,ione)
        if(matmodel.eq.1) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_1,td_matinit_1,prestr_mat_1,matchg,      ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.2) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_2,td_matinit_2,prestr_mat_2,matchg,      ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.3) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_3,td_matinit_3,prestr_mat_3,matchg,      ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.4) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_4,td_matinit_4,prestr_mat_4,matchg,      ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.5) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_5,td_matinit_5,prestr_mat_5,matchg,      ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.6) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_6,td_matinit_6,prestr_mat_6,matchg,      ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.7) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_7,td_matinit_7,prestr_mat_7,matchg,      ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.8) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_8,td_matinit_8,prestr_mat_8,matchg,      ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.9) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_9,td_matinit_9,prestr_mat_9,matchg,      ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.10) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_10,td_matinit_10,prestr_mat_10,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.11) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_11,td_matinit_11,prestr_mat_11,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.12) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_12,td_matinit_12,prestr_mat_12,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.13) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_13,td_matinit_13,prestr_mat_13,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.14) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_14,td_matinit_14,prestr_mat_14,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.15) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_15,td_matinit_15,prestr_mat_15,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.16) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_16,td_matinit_16,prestr_mat_16,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.17) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_17,td_matinit_17,prestr_mat_17,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.18) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_18,td_matinit_18,prestr_mat_18,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.19) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_19,td_matinit_19,prestr_mat_19,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.20) then
          call matinit_cmp(
     &     A,nnz,neq,                                                   ! sparse
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg), ! elemfamily
     &     nelfamily,nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states, ! elemfamily
     &     ielg,numelv,ielindx,                                         ! elemnt,orig element num
     &     ptmp,nprop,elas_mat_20,td_matinit_20,prestr_mat_20,matchg,   ! materl
     &     tminmax,                                                     ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else
          ierr=101
          errstrng="matinit_drv"
        end if
        if(ierr.ne.izero) return
      end do
c
c...  add Winkler elements to stiffness matrix diagonals
c
      if(nwink.ne.izero) call winklr(A,iwink,wink,histry,nstep,
     & nwink,nhist,nnz,lastep,ierr,errstrng)
      if(ierr.ne.izero) return
      if(nwinkx.ne.izero) call winklr(A,iwinkx,winkx,histry,nstep,
     & nwinkx,nhist,nnz,lastep,ierr,errstrng)
      if(ierr.ne.izero) return
c
c...  check stiffness matrix for zero or negative diagonals
c
      call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
      call ckdiag(A,ierr,errstrng)
      return
      end
c
c version
c $Id: matinit_drv.F,v 1.6 2005/05/05 21:15:00 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
