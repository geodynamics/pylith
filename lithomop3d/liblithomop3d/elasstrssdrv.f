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
      subroutine elstrssdrv(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     & numslp,
     & state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     & ien,lm,lmx,lmf,infiel,numelt,nconsz,
     & infmat,numat,
     & infetype,gauss,sh,shj,netypes,ngaussmax,nenmax,
     & skew,nskdim,numrot,
     & idebug,idout,kto,kw,getshape,bmatrix,ierr)
c
c...program to compute the total stress and strain for the current
c   iteration (elastic case).  This driver is for the small strain case.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nsd,ndof,numnp,neq,numfn,numslp,nstr,nddmat,nstatesz
      integer ndmatsz,numelt,nconsz,numat,netypes,ngaussmax,nenmax
      integer nskdim,numrot,idebug,idout,kto,kw,ierr
      integer ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz),lmf(nconsz)
      integer infiel(6,numelt),infmat(3,numat),infetype(4,netypes)
      double precision x(nsd,numnp),b(neq),d(ndof,numnp),dx(ndof,numnp)
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
      include "elasstrs_ext.inc"
      external getshape,bmatrix
c
c...  local variables
c
      integer matgpt,imat,matmodel,nmatel
c
cdebug      write(6,*) "Hello from elstrsss_f!"
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
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs1,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.2) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs2,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.3) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs3,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.4) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs4,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.5) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs5,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.6) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs6,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.7) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs7,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.8) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs8,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.9) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs9,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.10) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs10,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.11) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs11,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.12) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs12,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.13) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs13,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.14) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs14,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.15) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs15,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.16) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs16,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.17) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs17,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.18) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs18,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.19) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs19,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else if(matmodel.eq.20) then
          call elasstrsscmp(x,b,d,dx,tfault,nsd,ndof,numnp,neq,numfn,
     &     numslp,
     &     state,dmat,nstr,nddmat,nstatesz,ndmatsz,
     &     ien,lm,lmx,lmf,infiel,numelt,nconsz,
     &     infmat(1,imat),matgpt,elasstrs20,getshape,bmatrix,
     &     infetype,gauss,netypes,
     &     skew,nskdim,numrot,
     &     idebug,idout,kto,kw,ierr)
        else
          write(kto,*) "Error!  Material model does not exist!"
          stop
        end if
        if(ierr.ne.0) return
        matgpt=matgpt+nmatel
      end do
      return
      end
c
c version
c $Id: elasstrssdrv.f,v 1.1 2004/05/24 21:01:45 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
