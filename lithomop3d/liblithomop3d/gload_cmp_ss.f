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
      subroutine gload_cmp_ss(
     & gvec2,grav,neq,                                                  ! force
     & x,d,numnp,                                                       ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! fault
     & ien,lm,lmx,lmf,infiel,numelt,nconsz,                             ! elemnt
     & dens,infmat,matgpt,nmatel,matchg,                                ! materl
     & gauss,shj,infetype,                                              ! eltype
     & rtimdat,ntimdat,                                                 ! timdat
     & skew,numrot,                                                     ! skew
     & ierr,errstrng)                                                   ! errcode
c
c...  computation routine to add body forces due to gravity for small
c     strain problems
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
      integer neq,numnp,numslp,numfn,numelt,nconsz,matgpt,nmatel,numrot
      integer ierr
      integer ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz),lmf(nconsz)
      integer infiel(6,numelt),infmat(6),infetype(4,netypes)
      character errstrng*(*)
      logical matchg
      double precision gvec2(neq),grav(ndof),x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp),tfault(ndof,numfn),dens
      double precision gauss(nsd+1,ngaussmax,netypes)
      double precision shj(nsd+1,nenmax,ngaussmax,netypes)
      double precision skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
c
c...  local variables
c
      integer ind,iel,indien,ietype,ngauss,nen,nee
      double precision p(60),xl(60)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from gload_cmp_ss_f!"
c
c
c...  loop over elements in a material group
c
      do ind=matgpt,matgpt+nmatel-1
        iel=infiel(4,ind)
        indien=infiel(1,iel)
        ietype=infiel(3,iel)
        ngauss=infetype(1,ietype)
        nen=infetype(2,ietype)
        nee=infetype(4,ietype)
c
c...  localize coordinates
c
        call lcoord(x,xl,ien(indien),nen,numnp)
        call fill(p,zero,nee)
c
c...  compute element body force
c
        call gravld(p,grav,xl,iel,nen,dens,gauss(1,1,ietype),
     &   shj(1,1,1,ietype),ngauss,ierr,errstrng)
c
c...  rotate body forces if necessary
c
        if(numrot.ne.0) call rpforc(p,skew,ien(indien),numnp,nen)
c
c...  add element forces to global vector (gvec2)
c
        call addfor(gvec2,p,lm(1,indien),lmx(1,indien),neq,nee)
      end do
      return
      end
c
c version
c $Id: gload_cmp_ss.f,v 1.1 2004/07/01 20:18:46 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
