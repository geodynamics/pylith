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
      subroutine formf_ss(
     & bdisp,neq,                                                       ! force
     & x,nsd,ndof,numnp,                                                ! global
     & s,stemp,                                                         ! stiff
     & dmat,ien,lm,lmx,infiel,iddmat,nstr,nddmat,ndmatsz,numelt,nconsz, ! elemnt
     & infmat,numat,                                                    ! materl
     & gauss,sh,shj,infetype,netypes,                                   ! eltype
     & skew,nskdim,numrot,                                              ! skew
     & nfault,dfault,tfault,numfn,                                      ! split
     & getshape,bmatrix,                                                ! bbar
     & ierr)                                                            ! errcode
c
c      generates forces due to faulted nodes
c
      include "implicit.inc"
c
c...  dimension parameters
c
      include "nshape.inc"
c
c...  subroutine arguments
c
      integer neq,nsd,ndof,numnp,nstr,nddmat,ndmatsz,numelt,nconsz,numat
      integer netypes,nskdim,numrot,numfn,ierr
      integer ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz)
      integer infiel(6,numelt),iddmat(nstr,nstr),infmat(6,numat)
      integer infetype(4,netypes),nfault(3,numfn)
      double precision bdisp(neq),x(nsd,numnp),s(neemax*neemax)
      double precision stemp(neemax*neemax),dfault(ndof,numfn)
      double precision tfault(ndof,numfn),dmat(nddmat,ndmatsz)
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
c...  external routines
c
      external getshape,bmatrix
c
c...  local variables
c
      integer i,iel,indien,imat,imatvar,ietype,nen,inddmat,ngauss,nee
      integer ngaussdim
      double precision p(60),dl(60)
c
cdebug      write(6,*) "Hello from formf_ss_f!"
c
c
c...  loop over split node entries
c
      do i=1,numfn
        iel=nfault(1,i)
        indien=infiel(1,iel)
        imat=infiel(2,iel)
        imatvar=infmat(4,imat)
        ietype=infiel(3,iel)
        nen=infetype(2,ietype)
        inddmat=infiel(6,iel)
        ngauss=infetype(1,ietype)
        nee=infetype(4,ietype)
        ngaussdim=ngauss
        if(imatvar.eq.izero) ngaussdim=ngaussmax
c
c...  form element stiffness matrix
c
        call formes_ss(
     &   x,nsd,ndof,numnp,                                              ! global
     &   s,stemp,                                                       ! stiff
     &   dmat(1,inddmat),ien(indien),lm(1,indien),iddmat,nstr,nddmat,   ! elemnt
     &   iel,                                                           ! elemnt
     &   gauss(1,1,ietype),sh(1,1,1,ietype),shj(1,1,1,ietype),          ! eltype
     &   ngauss,ngaussdim,nen,nee,                                      ! eltype
     &   skew,nskdim,numrot,                                            ! skew
     &   getshape,bmatrix,                                              ! bbar
     &   ierr)                                                          ! errcode
c
        if(ierr.ne.izero) return
c
        call lflteq(dl,dfault(1,i),nfault(1,i),ien(indien),nen,ndof)
	call dsymv("u",nee,one,s,nee,dl,ione,zero,p,ione)
        call addfor(bdisp,p,lm(1,indien),lmx(1,indien),neq,nee)
      end do
      return
      end
c
c version
c $Id: formf_ss.f,v 1.1 2004/06/16 21:14:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
