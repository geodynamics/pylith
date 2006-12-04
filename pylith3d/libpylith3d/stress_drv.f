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
      subroutine stress_drv(
     & bintern,neq,                                                     ! force
     & x,d,numnp,iddmat,                                                ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! fault
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,ivfamily,nvfamilies,     ! elemnt
     & numelv,nstatesz,nstatesz0,nprestrflag,ipstrs,ipauto,             ! elemnt
     & prop,infmatmod,npropsz,tminmax,                                  ! materl
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
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
      include "materials.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer neq,numnp,numslp,numfn,nvfamilies,numelv,nstatesz
      integer nstatesz0,nprestrflag,ipstrs,ipauto
      integer npropsz,nen,ngauss,nee,nhist,lastep,numrot,ierr
      integer iddmat(nstr,nstr)
      integer ien(nen,numelv),lm(ndof*nen,numelv),lmx(ndof*nen,numelv)
      integer lmf(nen,numelv),ivfamily(6,nvfamilies)
      integer infmatmod(6,nmatmodmax)
      character errstrng*(*)
      double precision bintern(neq),x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp)
      double precision tfault(ndof,numfn)
      double precision state(nstatesz),dstate(nstatesz)
      double precision state0(nstatesz0)
      double precision dmat(nddmat*ngauss,numelv),prop(npropsz),tminmax
      double precision gauss(nsd+1,ngauss)
      double precision sh(nsd+1,nen,ngauss),shj(nsd+1,nen,ngauss)
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
      integer ielg,ifam,nelfamily,matmodel,indstate,indstate0,indprop
      integer nstate,nprop,nstate0,n0states
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
      tminmax=big
      ielg=ione
c
c...  loop over material groups and then select appropriate material model
c     routine
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
c*****************************
c       Temporary kludge.  In the near future, properties will be
c       obtained from a spatial database.  For now, an element family
c       is assumed to be defined by its material type.  In the next
c       refinement, it will be defined by a material model (and
c       possibly other characteristics such as element type), and
c       actual properties will be defined by the database for a given
c       location/time.
c*****************************
        matchg=.false.
        call dcopy(nprop,prop(indprop),ione,ptmp,ione)
        if(matmodel.eq.1) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_1,td_strs_1,matchg,tminmax,             ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.2) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_2,td_strs_2,matchg,tminmax,             ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.3) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_3,td_strs_3,matchg,tminmax,             ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.4) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_4,td_strs_4,matchg,tminmax,             ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.5) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_5,td_strs_5,matchg,tminmax,             ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.6) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_6,td_strs_6,matchg,tminmax,             ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.7) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_7,td_strs_7,matchg,tminmax,             ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.8) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_8,td_strs_8,matchg,tminmax,             ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.9) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_9,td_strs_9,matchg,tminmax,             ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.10) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_10,td_strs_10,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.11) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_11,td_strs_11,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.12) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_12,td_strs_12,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.13) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_13,td_strs_13,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.14) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_14,td_strs_14,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.15) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_15,td_strs_15,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.16) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_16,td_strs_16,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.17) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_17,td_strs_17,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.18) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_18,td_strs_18,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.19) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_19,td_strs_19,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.20) then
          call stress_cmp(
     &     bintern,neq,                                                 ! force
     &     x,d,numnp,iddmat,                                            ! global
     &     dx,numslp,                                                   ! slip
     &     tfault,numfn,                                                ! fault
     &     state(indstate),dstate(indstate),state0(indstate0),          ! elemfamily
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),lmx(1,ielg),             ! elemfamily
     &     lmf(1,ielg),nelfamily,nstate,nstate0,nprestrflag,ipstrs,     ! elemfamily
     &     ipauto,n0states,ielg,                                        ! elemfamily
     &     ptmp,nprop,elas_strs_20,td_strs_20,matchg,tminmax,           ! materl
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
     &     rtimdat,ntimdat,rgiter,                                      ! timdat
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
        else
          ierr=101
          errstrng="stress_drv"
        end if
        if(ierr.ne.izero) return
      end do
      return
      end
c
c version
c $Id: stress_drv.f,v 1.15 2005/04/16 00:40:50 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
