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
     & dmat,infmat,nstr,nddmat,ndmatsz,numat,                           ! materl
     & infiel,numelt,                                                   ! elmnt
     & prop,infprop,mhist,iddmat,npropsz,                               ! props
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
      integer nstr,nddmat,ndmatsz,numat,numelt,npropsz,netypes,ngaussmax
      integer nhist,nstep,lastep,idout,kto,kw
      integer infmat(4,numat),infiel(6,numelt),infprop(2,numat)
      integer mhist(2,numat),iddmat(nstr,nstr),infetype(4,netypes)
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
        nprop=infprop(1,imat)
        indprop=infprop(2,imat)
        matchg=.false.
        call mathist(ptmp,prop(indprop),mhist(1,indprop),histry,nprop,
     &   imat,nstep,nhist,lastep,matchg,idout,kto,kw)
        if(matmodel.eq.1) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat1,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.2) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat2,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.3) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat3,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.4) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat4,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.5) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat5,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.6) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat6,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.7) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat7,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.8) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat8,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.9) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat9,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.10) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat10,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.11) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat11,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.12) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat12,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.13) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat13,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.14) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat14,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.15) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat15,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.16) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat16,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.17) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat17,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.18) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat18,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.19) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat19,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        else if(matmodel.eq.20) then
          call elasmatcmp(dmat,ptmp,nstr,nddmat,ndmatsz,nprop,
     &     infmat(1,imat),matgpt,elasmat20,
     &     infiel,numelt,
     &     infetype,netypes,ngaussmax,
     &     iddmat)
        end if
        matgpt=matgpt+nmatel
      end do
      return
      end
c
c version
c $Id: elasmatdrv.f,v 1.1 2004/05/24 21:01:45 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
