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
      subroutine matinit_drv(
     & alnz,ja,nnz,neq,                                                 ! sparse
     & x,d,iwink,wink,numnp,nwink,                                      ! global
     & dx,iwinkx,winkx,numslp,numsn,nwinkx,                             ! slip
     & tfault,numfn,                                                    ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,ndmatsz, ! elemnt
     & numelt,nconsz,                                                   ! elemnt
     & prop,mhist,infmat,infmatmod,numat,npropsz,tminmax,               ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
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
      integer nnz,neq,numnp,nwink,numslp,numsn,nwinkx,numfn,nstatesz
      integer ndmatsz,numelt,nconsz,numat,npropsz,nhist,lastep,numrot
      integer ierr
      integer ja(nnz),iwink(2,nwink),iwinkx(2,nwinkx),ien(nconsz)
      integer lm(ndof,nconsz),lmx(ndof,nconsz),lmf(nconsz)
      integer infiel(6,numelt),iddmat(nstr,nstr),mhist(npropsz)
      integer infmat(3,numat),infmatmod(5,nmatmodmax)
      integer infetype(4,netypes)
      character errstrng*(*)
      double precision alnz(nnz),x(nsd,numnp),d(ndof,numnp)
      double precision wink(nwink),dx(ndof,numnp),winkx(nwinkx)
      double precision tfault(ndof,numfn)
      double precision s(neemax*neemax),stemp(neemax*neemax)
      double precision state(nstr,nstatesz),dstate(nstr,nstatesz)
      double precision dmat(nddmat,ndmatsz),prop(npropsz),tminmax
      double precision gauss(nsd+1,ngaussmax,netypes)
      double precision sh(nsd+1,nenmax,ngaussmax,netypes)
      double precision shj(nsd+1,nenmax,ngaussmax,netypes)
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
      external matinit_cmp,getshape,bmatrix
c
c...  local variables
c
      integer imat,matmodel,nmatel,imatvar,nstate,nprop,indprop,matgpt
      double precision ptmp(100)
      logical matchg
cdebug      integer idb,jdb
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
      include "rgiter_def.inc"
c
cdebug      write(6,*) "Hello from matinit_drv_f!"
c
      call fill(alnz,zero,nnz)
c*      call fill(dmat,zero,nddmat*ndmatsz)
      matgpt=1
      tminmax=big
c
c...  loop over material groups and then select appropriate material
c     model routine
c
      do imat=1,numat
        matmodel=infmat(1,imat)
        nmatel=infmat(2,imat)
        indprop=infmat(3,imat)
        nstate=infmatmod(2,matmodel)
        nprop=infmatmod(3,matmodel)
        imatvar=infmatmod(4,matmodel)
        matchg=.false.
        call mathist(ptmp,prop(indprop),mhist(indprop),histry,nprop,
     &   imat,nstep,nhist,lastep,matchg,ierr,errstrng)
        if(ierr.ne.izero) return
        if(matmodel.eq.1) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_1,          ! materl
     &     td_matinit_1,matchg,tminmax,                                 ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.2) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_2,          ! materl
     &     td_matinit_2,matchg,tminmax,                                 ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.3) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_3,          ! materl
     &     td_matinit_3,matchg,tminmax,                                 ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.4) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_4,          ! materl
     &     td_matinit_4,matchg,tminmax,                                 ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.5) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_5,          ! materl
     &     td_matinit_5,matchg,tminmax,                                 ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.6) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_6,          ! materl
     &     td_matinit_6,matchg,tminmax,                                 ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.7) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_7,          ! materl
     &     td_matinit_7,matchg,tminmax,                                 ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.8) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_8,          ! materl
     &     td_matinit_8,matchg,tminmax,                                 ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.9) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_9,          ! materl
     &     td_matinit_9,matchg,tminmax,                                 ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.10) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_10,         ! materl
     &     td_matinit_10,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.11) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_11,         ! materl
     &     td_matinit_11,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.12) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_12,         ! materl
     &     td_matinit_12,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.13) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_13,         ! materl
     &     td_matinit_13,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.14) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_14,         ! materl
     &     td_matinit_14,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.15) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_15,         ! materl
     &     td_matinit_15,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.16) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_16,         ! materl
     &     td_matinit_16,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.17) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_17,         ! materl
     &     td_matinit_17,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.18) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_18,         ! materl
     &     td_matinit_18,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.19) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_19,         ! materl
     &     td_matinit_19,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.20) then
          call matinit_cmp(
     &     alnz,ja,nnz,neq,                                             ! sparse
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,     ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     ptmp,nstate,nprop,matgpt,imatvar,nmatel,elas_mat_20,         ! materl
     &     td_matinit_20,matchg,tminmax,                                ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else
          ierr=101
          errstrng="matinit_drv"
        end if
        if(ierr.ne.izero) return
        matgpt=matgpt+nmatel
      end do
c
c...  debugging output
c
cdebug      write(6,*) "DMAT array:"
cdebug      do idb=1,ndmatsz
cdebug        write(6,*) (dmat(jdb,idb),jdb=1,nddmat)
cdebug      end do
cdebug      write(6,*) "ALNZ array:"
cdebug      write(6,*) (alnz(idb),idb=1,nnz)
c
c...  add Winkler elements to stiffness matrix diagonals
c
      if(nwink.ne.izero) call winklr(alnz,iwink,wink,histry,nstep,
     & nwink,nhist,nnz,lastep,ierr,errstrng)
      if(ierr.ne.izero) return
      if(nwinkx.ne.izero) call winklr(alnz,iwinkx,winkx,histry,nstep,
     & nwinkx,nhist,nnz,lastep,ierr,errstrng)
      if(ierr.ne.izero) return
c
c...  check stiffness matrix for zero or negative diagonals
c
      call ckdiag(alnz,neq,nnz,ierr,errstrng)
      return
      end
c
c version
c $Id: matinit_drv.f,v 1.5 2004/08/12 02:07:33 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
