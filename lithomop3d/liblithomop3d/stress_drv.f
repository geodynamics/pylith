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
      subroutine stress_drv(
     & b,neq,                                                           ! force
     & x,d,numnp,                                                       ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! fault
     & state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,numelt, ! elemnt
     & nconsz,                                                          ! elemnt
     & prop,mhist,infmat,numat,npropsz,                                 ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & histry,rtimdat,rgiter,ntimdat,nhist,lastep,stress_cmp,           ! timdat
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...program to compute the total stress and strain for the current
c   iteration and compute the internal force vector.
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
      integer neq,numnp,numslp,numfn,nstatesz,ndmatsz,numelt,nconsz
      integer numat,npropsz,nhist,lastep,numrot,ierr
      integer ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz),lmf(nconsz)
      integer infiel(6,numelt),mhist(npropsz),infmat(6,numat)
      integer infetype(4,netypes)
      character errstrng*(*)
      double precision b(neq),x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn)
      double precision state(nstr,nstatesz),dstate(nstr,nstatesz)
      double precision dmat(nddmat,ndmatsz),prop(npropsz)
      double precision gauss(nsd+1,ngaussmax,netypes)
      double precision sh(nsd+1,nenmax,ngaussmax,netypes)
      double precision shj(nsd+1,nenmax,ngaussmax,netypes)
      double precision histry(nhist,lastep+1),skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "ntimdat_dim.inc"
c
c...  intrinsic functions
c
c
c...  external routines
c
      include "elas_strs_ext.inc"
      include "td_strs_ext.inc"
      external stress_cmp,getshape,bmatrix
c
c...  local variables
c
      integer matgpt,imat,matmodel,nmatel,nprop,indprop
      logical matchg
      double precision ptmp(100)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from stress_drv_f!"
c
      call fill(b,zero,neq)
      matgpt=1
c
c...  loop over material groups and then select appropriate material model
c     routine
c
      do imat=1,numat
        matmodel=infmat(1,imat)
        nmatel=infmat(2,imat)
        nprop=infmat(5,imat)
        indprop=infmat(6,imat)
        matchg=.false.
        call mathist(ptmp,prop(indprop),mhist(indprop),histry,nprop,
     &   imat,nstep,nhist,lastep,matchg,ierr,errstrng)
        if(ierr.ne.izero) return
        if(matmodel.eq.1) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_1,td_strs_1,      ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.2) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_2,td_strs_2,      ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.3) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_3,td_strs_3,      ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.4) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_4,td_strs_4,      ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.5) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_5,td_strs_5,      ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.6) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_6,td_strs_6,      ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.7) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_7,td_strs_7,      ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.8) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_8,td_strs_8,      ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.9) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_9,td_strs_9,      ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.10) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_10,td_strs_10,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.11) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_11,td_strs_11,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.12) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_12,td_strs_12,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.13) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_13,td_strs_13,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.14) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_14,td_strs_14,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.15) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_15,td_strs_15,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.16) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_16,td_strs_16,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.17) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_17,td_strs_17,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.18) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_18,td_strs_18,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.19) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_19,td_strs_19,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.20) then
          call stress_cmp(
     &     b,neq,                                                       ! force
     &     x,d,numnp,                                                   ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,    ! elemnt
     &     numelt,nconsz,                                               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_strs_20,td_strs_20,    ! materl
     &     gauss,sh,shj,infetype,                                       ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else
          ierr=101
          errstrng="stress_drv"
        end if
        if(ierr.ne.izero) return
        matgpt=matgpt+nmatel
      end do
      return
      end
c
c version
c $Id: stress_drv.f,v 1.3 2004/07/02 18:07:28 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
