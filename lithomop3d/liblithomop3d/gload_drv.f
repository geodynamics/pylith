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
      subroutine gload_drv(
     & bgravity,ngravflag,grav,neq,                                     ! force
     & x,d,numnp,                                                       ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! fault
     & ien,lm,lmx,lmf,ivfamily,nvfamilies,numelv,                       ! elemnt
     & prop,mhist,infmatmod,npropsz,                                    ! materl
     & gauss,shj,nen,ngauss,nee,                                        ! eltype
     & histry,rtimdat,ntimdat,nhist,lastep,gload_cmp,                   ! timdat
     & skew,numrot,                                                     ! skew
     & ierr,errstrng)                                                   ! errcode
c
c...  driver routine to add body forces due to gravity
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer ngravflag,neq,numnp,numslp,numfn,nvfamilies,numelv,npropsz
      integer nen,ngauss,nee,nhist,lastep,numrot,ierr
      integer ien(nen,numelv),lm(ndof*nen,numelv),lmx(ndof*nen,numelv)
      integer lmf(nen,numelv),ivfamily(5,nvfamilies)
      integer mhist(npropsz),infmatmod(6,nmatmodmax)
      character errstrng*(*)
      double precision bgravity(ngravflag*neq),grav(ndof)
      double precision x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn),prop(npropsz)
      double precision gauss(nsd+1,ngauss)
      double precision shj(nsd+1,nen,ngauss)
      double precision histry(nhist,lastep+1),skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
c
c...  external routines
c
      external gload_cmp
c
c...  local variables
c
      integer i,ielg,ifam,nelfamily,matmodel,indprop,nprop,imat
      logical matchg
      double precision gsum,dens
      double precision ptmp(100)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from gload_drv_f!"
c
      gsum=zero
      do i=1,ndof
        gsum=gsum+grav(i)*grav(i)
      end do
      if(gsum.eq.zero) return
      call fill(bgravity,zero,neq)
      ielg=ione
c
c...  loop over element families
c
      do ifam=1,nvfamilies
        nelfamily=ivfamily(1,ifam)
        matmodel=ivfamily(2,ifam)
        indprop=ivfamily(5,ifam)
        nprop=infmatmod(3,matmodel)
        imat=ifam
        matchg=.false.
        call mathist(ptmp,prop(indprop),mhist(indprop),histry,nprop,
     &   imat,nstep,nhist,lastep,matchg,ierr,errstrng)
c
        if(ierr.ne.izero) return
c
        dens=ptmp(1)
c
        call gload_cmp(
     &   bgravity,grav,neq,                                             ! force
     &   x,d,numnp,                                                     ! global
     &   dx,numslp,                                                     ! slip
     &   tfault,numfn,                                                  ! fault
     &   ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg),nelfamily,ielg, ! elemfamily
     &   dens,matchg,                                                   ! materl
     &   gauss,shj,nen,ngauss,nee,                                      ! eltype
     &   rtimdat,ntimdat,                                               ! timdat
     &   skew,numrot,                                                   ! skew
     &   ierr,errstrng)                                                 ! errcode
c
      end do
c
      return
      end
c
c version
c $Id: gload_drv.f,v 1.5 2005/03/21 22:13:15 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
