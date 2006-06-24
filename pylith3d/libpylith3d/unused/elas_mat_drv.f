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
      subroutine elas_mat_drv(
     & dmat,infiel,iddmat,ndmatsz,numelt,                               ! elemnt
     & prop,mhist,infmat,numat,npropsz,                                 ! materl
     & infetype,                                                        ! eltype
     & histry,nhist,nstep,lastep,                                       ! timdat
     & ierr,errstrng)                                                   ! errcode
c
c...  driver subroutine to form the d-matrix for the elastic solution
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
      integer ndmatsz,numelt,numat,npropsz
      integer nhist,nstep,lastep,ierr
      integer infiel(6,numelt),iddmat(nstr,nstr),mhist(npropsz)
      integer infmat(6,numat),infetype(4,netypes)
      character errstrng*(*)
      double precision dmat(nddmat,ndmatsz),prop(npropsz)
      double precision histry(nhist,lastep+1)
c
c...  external routines
c
      include "elas_mat_ext.inc"
c
c...  local variables
c
      integer imat,matmodel,nmatel,nprop,indprop,matgpt
      double precision ptmp(100)
      logical matchg
c
cdebug      write(6,*) "Hello from elas_mat_drv_f!"
c
      call fill(dmat,zero,nddmat*ndmatsz)
      matgpt=1
c
c...  loop over material groups and then select appropriate material
c     model routine
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
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_1,                 ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.2) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_2,                 ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.3) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_3,                 ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.4) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_4,                 ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.5) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_5,                 ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.6) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_6,                 ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.7) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_7,                 ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.8) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_8,                 ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.9) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_9,                 ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.10) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_10,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.11) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_11,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.12) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_12,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.13) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_13,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.14) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_14,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.15) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_15,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.16) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_16,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.17) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_17,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.18) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_18,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.19) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_19,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else if(matmodel.eq.20) then
          call elas_mat_cmp(
     &     dmat,infiel,iddmat,ndmatsz,numelt,                           ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elas_mat_20,                ! materl
     &     infetype,                                                    ! eltype
     &     ierr,errstrng)                                               ! errcode
        else
          ierr=101
          errstrng="elas_mat_drv"
        end if
        if(ierr.ne.izero) return
        matgpt=matgpt+nmatel
      end do
      return
      end
c
c version
c $Id: elas_mat_drv.f,v 1.1 2004/06/25 21:41:22 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
