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
      subroutine formes_ss(
     & x,numnp,iddmat,                                                  ! global
     & s,stemp,                                                         ! stiff
     & dmat,ien,lm,iel,                                                 ! elemnt
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...  subroutine to form the elemental stiffness matrix
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
      integer numnp,iel,nen,ngauss,nee,numrot,ierr
      integer iddmat(nstr,nstr),ien(nen),lm(ndof,nen)
      character errstrng*(*)
      double precision x(nsd,numnp),s(neemax*neemax)
      double precision stemp(neemax*neemax),dmat(nddmat,ngauss)
      double precision gauss(nsd+1,ngauss),sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss),skew(nskdim,numnp)
c
c...  external routines
c
      external getshape,bmatrix
c
c...  local variables
c
      double precision xl(60)
c
cdebug      write(6,*) "Hello from formes_ss_f!"
c
      call fill(s,zero,nee*nee)
c
c...  localize coordinates
c
      call lcoord(x,xl,ien,nen,numnp)
c
c...  construct local stiffness matrix, symmetrize it, and rotate for
c     skew boundary conditions
c
      call stiff_ss(
     & xl,iddmat,                                                       ! global
     & dmat,ien,iel,                                                    ! elemnt
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & s,                                                               ! stiff
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
      if(ierr.ne.izero) return
c
      if(numrot.ne.izero) call rstiff(s,stemp,skew,ien,numnp,nen,nee)
c
      return
      end
c
c version
c $Id: formes_ss.f,v 1.6 2005/04/01 23:34:13 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
