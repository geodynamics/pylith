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
      subroutine elasmatdrv(
     & dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,                   ! elemnt
     & prop,mhist,infmat,numat,npropsz,                                 ! materl
     & infetype,netypes,ngaussmax,                                      ! eltype
     & histry,nhist,nstep,lastep,                                       ! timdat
     & idout,kto,kw)                                                    ! fileinf
c
c...  driver subroutine to form the d-matrix for the elastic solution
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,nddmat,ndmatsz,numelt,numat,npropsz,netypes,ngaussmax
      integer nhist,nstep,lastep,idout,kto,kw
      integer infiel(6,numelt),iddmat(nstr,nstr),mhist(npropsz)
      integer infmat(6,numat),infetype(4,netypes)
      double precision dmat(nddmat,ndmatsz),prop(npropsz)
      double precision histry(nhist,lastep+1)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  external routines
c
      include "elasmat_ext.inc"
c
c...  local variables
c
      integer imat,matmodel,nmatel,nprop,indprop,matgpt
      double precision ptmp(100)
      logical matchg
c
cdebug      write(6,*) "Hello from elasmatdrv_f!"
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
     &   imat,nstep,nhist,lastep,matchg,idout,kto,kw)
        if(matmodel.eq.1) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat1,                   ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.2) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat2,                   ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.3) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat3,                   ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.4) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat4,                   ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.5) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat5,                   ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.6) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat6,                   ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.7) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat7,                   ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.8) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat8,                   ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.9) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat9,                   ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.10) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat10,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.11) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat11,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.12) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat12,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.13) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat13,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.14) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat14,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.15) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat15,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.16) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat16,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.17) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat17,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.18) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat18,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.19) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat19,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        else if(matmodel.eq.20) then
          call elasmatcmp(
     &     dmat,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,               ! elemnt
     &     ptmp,infmat(1,imat),nprop,matgpt,elasmat20,                  ! materl
     &     infetype,netypes,ngaussmax)                                  ! eltype
        end if
        matgpt=matgpt+nmatel
      end do
      return
      end
c
c version
c $Id: elasmatdrv.f,v 1.3 2004/05/25 17:38:55 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
