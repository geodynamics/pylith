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
      subroutine gload_cmp_ss(
     & bgravity,grav,neq,                                               ! force
     & x,d,numnp,                                                       ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! fault
     & ien,lm,lmx,lmf,nelfamily,ielg,                                   ! elemfamily
     & dens,matchg,                                                     ! materl
     & gauss,shj,nen,ngauss,nee,                                        ! eltype
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
      integer neq,numnp,numslp,numfn,nelfamily,ielg,nen,ngauss,nee
      integer numrot,ierr
      integer ien(nen,nelfamily),lm(ndof*nen,nelfamily)
      integer lmx(ndof*nen,nelfamily),lmf(nen,nelfamily)
      character errstrng*(*)
      logical matchg
      double precision bgravity(neq),grav(ndof)
      double precision x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp),tfault(ndof,numfn),dens
      double precision gauss(nsd+1,ngauss)
      double precision shj(nsd+1,nen,ngauss)
      double precision skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
c
c...  local variables
c
      integer ielf
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
c...  loop over elements in a family
c
      do ielf=1,nelfamily
c
c...  localize coordinates
c
        call lcoord(x,xl,ien(1,ielf),nen,numnp)
        call fill(p,zero,nee)
c
c...  compute element body force
c
        call gravld(p,grav,xl,ielg,nen,dens,gauss,shj,ngauss,ierr,
     &   errstrng)
c
c...  rotate body forces if necessary
c
        if(numrot.ne.izero) call rpforc(p,skew,ien(1,ielf),numnp,nen)
c
c...  add element forces to global vector (bgravity)
c
        call addfor(bgravity,p,lm(1,ielf),lmx(1,ielf),neq,nee)
        ielg=ielg+ione
      end do
      return
      end
c
c version
c $Id: gload_cmp_ss.f,v 1.6 2005/04/08 00:41:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
