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
      subroutine stress_mat_drv(
     & A,nnz,                                                           ! sparse
     & bintern,neq,                                                     ! force
     & x,d,iwink,wink,numnp,nwink,iddmat,                               ! global
     & dx,iwinkx,winkx,numslp,numsn,nwinkx,                             ! slip
     & tfault,numfn,                                                    ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,nvfamilies,     ! elemnt
     & numelv,nstatesz,nstatesz0,nprestrflag,ipstrs,ipauto,             ! elemnt
     & prop,mhist,infmatmod,npropsz,tminmax,                            ! materl
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & histry,rtimdat,rgiter,ntimdat,nhist,lastep,stress_mat_cmp,       ! timdat
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...program to compute the total stress and strain for the current
c   iteration and compute the internal force vector.
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
      integer iddmat(nstr,nstr)
      integer iwink(2,nwink),iwinkx(2,nwinkx),ien(nen,numelv)
      integer lm(ndof*nen,numelv),lmx(ndof*nen,numelv),lmf(nen,numelv)
      integer ivfamily(5,nvfamilies),mhist(npropsz)
      integer infmatmod(6,nmatmodmax)
      character errstrng*(*)
      double precision bintern(neq),x(nsd,numnp),d(ndof,numnp)
      double precision wink(nwink),dx(ndof,numnp),winkx(nwinkx)
      double precision tfault(ndof,numfn),s(neemax*neemax)
      double precision stemp(neemax*neemax),state(nstatesz)
      double precision dstate(nstatesz),state0(nstatesz0)
      double precision dmat(nddmat*ngauss,numelv),prop(npropsz),tminmax
      double precision gauss(nsd+1,ngauss)
      double precision sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss)
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
      include "td_strs_mat_ext.inc"
      external stress_mat_cmp,getshape,bmatrix
c
c...  local variables
c
      integer ielg,ifam,nelfamily,matmodel,indstate,indstate0,indprop
      integer nstate,nprop,nstate0,n0states,imat
      logical matchg
      double precision ptmp(100)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from stress_mat_drv_f!"
c
      tminmax=big
      nrftot=ntimdat(7)
      nrftot=nrftot+1
      ntimdat(7)=nrftot
      ielg=ione
c
c...  loop over element families and then select appropriate material
c     model routine.
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
        call mathist(ptmp,prop(indprop),mhist(indprop),histry,nprop,
     &   imat,nstep,nhist,lastep,matchg,ierr,errstrng)
        if(ierr.ne.izero) return
c
        if(matmodel.eq.1) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_1,td_strs_mat_1,matchg,tminmax,         ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.2) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_2,td_strs_mat_2,matchg,tminmax,         ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.3) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_3,td_strs_mat_3,matchg,tminmax,         ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.4) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_4,td_strs_mat_4,matchg,tminmax,         ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.5) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_5,td_strs_mat_5,matchg,tminmax,         ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.6) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_6,td_strs_mat_6,matchg,tminmax,         ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.7) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_7,td_strs_mat_7,matchg,tminmax,         ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.8) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_8,td_strs_mat_8,matchg,tminmax,         ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.9) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_9,td_strs_mat_9,matchg,tminmax,         ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.10) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_10,td_strs_mat_10,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.11) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_11,td_strs_mat_11,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.12) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_12,td_strs_mat_12,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.13) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_13,td_strs_mat_13,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.14) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_14,td_strs_mat_14,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.15) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_15,td_strs_mat_15,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.16) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_16,td_strs_mat_16,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.17) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_17,td_strs_mat_17,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.18) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_18,td_strs_mat_18,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.19) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_19,td_strs_mat_19,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.20) then
          call stress_mat_cmp(
     &     A,nnz,                                                       ! sparse
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,numsn,                                             ! slip
     &     tfault,numfn,                                                ! fault
     &     s,stemp,                                                     ! stiff
     &     state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,           ! elemnt
     &     nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,      ! elemnt
     &     ptmp,nprop,elas_strs_20,td_strs_mat_20,matchg,tminmax,       ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else
          ierr=101
          errstrng="stress_mat_drv"
        end if
        if(ierr.ne.izero) return
      end do
      call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
      call ckdiag(A,ierr,errstrng)
      return
      end
c
c version
c $Id: stress_mat_drv.F,v 1.2 2005/03/22 19:46:09 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
