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
      subroutine elas_matinit_cmp_ss(
     & alnz,ja,nnz,neq,                                                 ! sparse
     & x,d,numnp,                                                       ! global
     & dx,numslp,numsn,                                                 ! slip
     & tfault,numfn,                                                    ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,         ! elemnt
     & ndmatsz,numelt,nconsz,                                           ! elemnt
     & prop,infmat,nprop,matgpt,elas_matinit,td_matinit,matchg,         ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & rtimdat,ntimdat,                                                 ! timdat
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...  compute subroutine to form the d-matrix for the elastic solution
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
      integer nnz,neq,numnp,numslp,numsn,numfn,nstatesz,ndmatsz,numelt
      integer nconsz,nprop,matgpt,numrot,ierr
      logical matchg
      integer ja(nnz),ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz)
      integer lmf(nconsz),infiel(6,numelt),iddmat(nstr,nstr),infmat(6)
      integer infetype(4,netypes)
      character errstrng*(*)
      double precision alnz(nnz),x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp),tfault(ndof,numfn)
      double precision s(neemax*neemax),stemp(neemax*neemax)
      double precision state(nstr,nstatesz),dstate(nstr,nstatesz)
      double precision dmat(nddmat,ndmatsz),prop(nprop)
      double precision gauss(nsd+1,ngaussmax,netypes)
      double precision sh(nsd+1,nenmax,ngaussmax,netypes)
      double precision shj(nsd+1,nenmax,ngaussmax,netypes)
      double precision skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
c
c...  external routines
c
      external elas_matinit,td_matinit,getshape,bmatrix
c
c...  local variables
c
      integer nmatel,imatvar,iel,inddmat0,ietype,ngauss,ng,inddmatg
      integer l,ind,inddmat,indien,nen,nee,ngtest,ngaussdim
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from elas_matinit_cmp_ss_f!"
c
      nmatel=infmat(2)
      imatvar=infmat(4)
c
c...  compute d-matrix for first element in the group.  If imatvar is
c     zero, there is only one d-matrix for the group, and that is all
c     that needs to be done.  Otherwise, since this is the elastic
c     solution, the same matrix can be copied into each slot.
c
      iel=infiel(4,matgpt)
      inddmat0=infiel(6,iel)
      ietype=infiel(3,iel)
      call elas_matinit(dmat(1,inddmat0),prop,iddmat,nprop,ierr,
     & errstrng)
      if(ierr.ne.izero) return
      ngauss=infetype(1,ietype)
      ng=ngauss
      if(imatvar.eq.izero) ng=ngaussmax
      inddmatg=inddmat0
      do l=2,ng
        inddmatg=inddmatg+nddmat
        call dcopy(nddmat,dmat(1,inddmat0),ione,dmat(1,inddmatg),ione)
      end do
c
c...  loop over elements in group if there is material property
c     variation for the material type.
c
      if(imatvar.ne.0) then
        do ind=matgpt+1,matgpt+nmatel-1
          iel=infiel(4,ind)
          ietype=infiel(3,iel)
          inddmat=infiel(6,iel)
          ngauss=infetype(1,ietype)
          inddmatg=inddmat
          do l=1,ngauss
            inddmatg=inddmat+nddmat
            call dcopy(nddmat,dmat(1,inddmat0),ione,dmat(1,inddmatg),
     &       ione)
          end do
        end do
      end if
c
c...  loop over elements in group and add element stiffness to global
c     stiffness.
c
      ngtest=0
      if(imatvar.eq.izero) ngtest=ngaussmax
      do ind=matgpt,matgpt+nmatel-1
        iel=infiel(4,ind)
        indien=infiel(1,iel)
        ietype=infiel(3,iel)
        nen=infetype(2,ietype)
        inddmat=infiel(6,iel)
        ngauss=infetype(1,ietype)
        nee=infetype(4,ietype)
        ngaussdim=max(ngauss,ngtest)
        call formes_ss(
     &   x,numnp,                                                       ! global
     &   s,stemp,                                                       ! stiff
     &   dmat(1,inddmat),ien(indien),lm(1,indien),iddmat,iel,           ! elemnt
     &   gauss(1,1,ietype),sh(1,1,1,ietype),shj(1,1,1,ietype),          ! eltype
     &   ngauss,ngaussdim,nen,nee,                                      ! eltype
     &   skew,numrot,                                                   ! skew
     &   getshape,bmatrix,                                              ! bbar
     &   ierr,errstrng)                                                 ! errcod
c
        if(ierr.ne.izero) return
c
        call addstf(alnz,s,lm(1,indien),lmx(1,indien),ja,nee,numsn,nnz)
      end do
c
      return
      end
c
c version
c $Id: elas_matinit_cmp_ss.f,v 1.1 2004/06/23 19:19:21 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
