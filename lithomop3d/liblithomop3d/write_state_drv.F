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
      subroutine write_state_drv(
     & state,dstate,ivfamily,nvfamilies,numelv,nstatesz,                ! elemnt
     & infmatmod,                                                       ! materl
     & ngauss,                                                          ! eltype
     & delt,nstep,                                                      ! timdat
     & istatout,nstatout,                                               ! ioopts
     & idout,idsk,iucd,kw,kp,kucd,ucdroot,iprestress,                   ! ioinfo
     & ierr,errstrng)                                                   ! errcode
c
c...  program to print state variables
c     Note:  at present, it is assumed that the same istatout array
c     is used for the elastic and time-dependent solutions.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "materials.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nvfamilies,numelv,nstatesz,ngauss,nstep
      integer idout,idsk,iucd,kw,kp,kucd,iprestress,ierr
      integer ivfamily(5,nvfamilies),infmatmod(6,nmatmodmax)
      integer istatout(nstatesmax,3),nstatout(3)
      double precision delt
      double precision state(nstatesz),dstate(nstatesz)
      character ucdroot*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  external routines
c
      include "getstate_ext.inc"
c
c...  local variables
c
      integer istatoutc(nstatesmax*3),ibyteg(3*nstatesmax)
      integer nstatestot,nout,npts,ind,i,j,iopt
      integer ifam,nelfamily,matmodel,indstate,nstate,ielg
      integer ibyte,intlen,floatlen,istride
      double precision delti
      double precision statemin(3*nstatesmax),statemax(3*nstatesmax)
c
c...  included variable definitions
c
c
cdebug      write(6,*) "Hello from write_state_drv_f!"
c
      if(delt.gt.zero) then
        delti=one/delt
      else
        delti=zero
      end if
      nstatestot=nstatout(1)+nstatout(2)+nstatout(3)
      if(nstatestot.eq.izero) return
      nout=izero
      npts=izero
c
c...  create compacted version of istatout array
c
      ind=izero
      do i=1,3
        do j=1,nstatout(i)
          ind=ind+ione
          istatoutc(ind)=istatout(j,i)+(i-1)*nstatesmax
cdebug          write(6,*) "i,j,ind,nstatout(i),istatoutc(ind),istatout(j,i):"
cdebug          write(6,*) i,j,ind,nstatout(i),istatoutc(ind),istatout(j,i)
        end do
      end do
c
c...  create and open UCD file if UCD output is desired
c
      if(iucd.ne.izero) then
        iopt=4
        call open_ucd(kucd,iprestress,nstep,ucdroot,iopt,iucd)
        call write_ucd_header(istatoutc,nstatestot,kucd,iucd)
      end if
c
c...  initialize max, min, and byte location values for binary UCD
c
      intlen=ifour
      floatlen=ifour
      call fill(statemin,big,3*nstatesmax)
      call fill(statemax,-big,3*nstatesmax)
      ibyte=ione+2048+intlen*(1+nstatestot)
      istride=ngauss*numelv*floatlen
      ibyteg(1)=ibyte+2*floatlen*nstatestot
      do i=2,nstatestot
        ibyteg(i)=ibyteg(i-1)+istride
      end do
c
c...  loop over element families
c
      ielg=ione
      do ifam=1,nvfamilies
        nelfamily=ivfamily(1,ifam)
        matmodel=ivfamily(2,ifam)
        indstate=ivfamily(3,ifam)
        nstate=infmatmod(2,matmodel)
c
        if(matmodel.eq.1) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_1,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.2) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_2,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.3) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_3,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.4) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_4,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.5) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_5,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.6) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_6,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.7) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_7,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.8) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_8,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.9) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_9,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.10) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_10,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.11) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_11,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.12) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_12,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.13) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_13,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.14) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_14,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.15) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_15,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.16) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_16,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.17) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_17,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.18) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_18,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.19) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_19,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else if(matmodel.eq.20) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,ielg,      ! elemfamily
     &     get_state_20,                                                ! materl
     &     ngauss,                                                      ! eltype
     &     delti,nstep,                                                 ! timdat
     &     istatout,istatoutc,nstatout,nstatestot,nout,npts,            ! ioopts
     &     statemin,statemax,ibyteg,                                    ! binucd
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
        else
          ierr=101
          errstrng="write_state_drv"
          return
        end if
        ielg=ielg+nelfamily
      end do
c
c...  output min and max values and footer for binary UCD file
c
      if(iucd.eq.2) then
c
c...  The current setup uses autoconf to determine whether the fortran
c     compiler supports old-style f77 unformatted direct access,
c     stream I/O, or neither.  If neither is supported, binary output
c     is turned off.
c
#if defined FORTRAN_F77_IO
        write(kucd,rec=ibyte) (real(statemin(istatoutc(i))),
     &   i=1,nstatestot)
        ibyte=ibyte+nstatestot*floatlen
        write(kucd,rec=ibyte) (real(statemax(istatoutc(i))),
     &   i=1,nstatestot)
        write(kucd,rec=ibyteg(nstatestot))
     &   (real(statemax(istatoutc(i))),i=1,nstatestot)
#elif defined FORTRAN_STREAM_IO
        write(kucd,pos=ibyte) (real(statemin(istatoutc(i))),
     &   i=1,nstatestot)
        ibyte=ibyte+nstatestot*floatlen
        write(kucd,pos=ibyte) (real(statemax(istatoutc(i))),
     &   i=1,nstatestot)
        write(kucd,pos=ibyteg(nstatestot))
     &   (real(statemax(istatoutc(i))),i=1,nstatestot)
#else
#warning "Binary UCD output not supported with gfortran!"
#endif
      end if
      if(iucd.gt.izero) close(kucd)
      return
      end
c
c version
c $Id: write_state_drv.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
