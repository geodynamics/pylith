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
      subroutine elas_strs_drv_ss(
     & b,neq,                                                           ! force
     & x,d,dx,nsd,ndof,numnp,                                           ! global
     & tfault,numfn,numslp,                                             ! fault
     & state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,ndmatsz,   ! elemnt
     & numelt,nconsz,                                                   ! elemnt
     & infmat,numat,                                                    ! materl
     & gauss,sh,shj,infetype,netypes,                                   ! eltype
     & skew,nskdim,numrot,                                              ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr)                                                            ! errcode
c
c...program to compute the total stress and strain for the current
c   iteration (elastic case).  This driver is for the small strain case.
c
      include "implicit.inc"
c
c...  dimension parameters
c
      include "nshape.inc"
c
c...  subroutine arguments
c
      integer nsd,ndof,numnp,neq,numfn,numslp,nstr,nddmat,nstatesz
      integer ndmatsz,numelt,nconsz,numat,netypes
      integer nskdim,numrot,ierr
      integer ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz),lmf(nconsz)
      integer infiel(6,numelt),infmat(6,numat),infetype(4,netypes)
      double precision b(neq),x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn)
      double precision state(nstr,nstatesz),dmat(nddmat,ndmatsz)
      double precision gauss(nsd+1,ngaussmax,netypes)
      double precision sh(nsd+1,nenmax,ngaussmax,netypes)
      double precision shj(nsd+1,nenmax,ngaussmax,netypes)
      double precision skew(nskdim,numnp)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  intrinsic functions
c
c
c...  external routines
c
      include "elas_strs_ext.inc"
      external getshape,bmatrix
c
c...  local variables
c
      integer matgpt,imat,matmodel,nmatel
c
cdebug      write(6,*) "Hello from elas_strs_drv_ss_f!"
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
        if(matmodel.eq.1) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_1,                           ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.2) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_2,                           ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.3) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_3,                           ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.4) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_4,                           ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.5) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_5,                           ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.6) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_6,                           ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.7) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_7,                           ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.8) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_8,                           ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.9) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_9,                           ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.10) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_10,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.11) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_11,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.12) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_12,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.13) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_13,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.14) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_14,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.15) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_15,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.16) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_16,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.17) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_17,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.18) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_18,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.19) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_19,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else if(matmodel.eq.20) then
          call elas_strs_cmp_ss(
     &     b,neq,                                                       ! force
     &     x,d,dx,nsd,ndof,numnp,                                       ! global
     &     tfault,numfn,numslp,                                         ! fault
     &     state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,       ! elemnt
     &     ndmatsz,numelt,nconsz,                                       ! elemnt
     &     infmat(1,imat),matgpt,elas_strs_20,                          ! materl
     &     gauss,sh,shj,infetype,netypes,                               ! eltype
     &     skew,nskdim,numrot,                                          ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr)                                                        ! errcode
        else
          ierr=101
        end if
        if(ierr.ne.0) return
        matgpt=matgpt+nmatel
      end do
      return
      end
c
c version
c $Id: elas_strs_drv_ss.f,v 1.2 2004/06/17 19:05:11 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
