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
        write(kucd,rec=ibyte) (real(statemin(istatoutc(i))),
     &   i=1,nstatestot)
        ibyte=ibyte+nstatestot*floatlen
        write(kucd,rec=ibyte) (real(statemax(istatoutc(i))),
     &   i=1,nstatestot)
        write(kucd,rec=ibyteg(nstatestot))
     &   (real(statemax(istatoutc(i))),i=1,nstatestot)
      end if
      if(iucd.gt.izero) close(kucd)
      return
      end
c
c version
c $Id: write_state_drv.f,v 1.1 2005/08/05 19:58:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
