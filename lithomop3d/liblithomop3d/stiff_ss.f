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
      subroutine stiff_ss(
     & xl,nsd,ndof,                                                     ! global
     & dmat,ien,nstr,nddmat,iel,                                        ! elemnt
     & gauss,sh,shj,ngauss,nen,ngaussdim,                               ! eltype
     & s,nee,                                                           ! stiff
     & getshape,bmatrix,                                                ! bbar
     & idout,kto,kw,                                                    ! ioinf
     & ierr)                                                            ! errcode
c
c...computes the local stiffness matrix at the given integration points.
c   k=(b)t*d*b.
c
      include "implicit.inc"
c
c...  dimension parameters
c
      include "nshape.inc"
c
c...  subroutine arguments
c
      integer nsd,ndof,nstr,nddmat,iel,ngauss,nen,nee,idout,kto,kw,ierr
      integer ien(nen)
      double precision xl(nsd,nen),dmat(nddmat,ngaussdim)
      double precision gauss(nsd+1,ngaussmax),sh(nsd+1,nenmax,ngaussmax)
      double precision shj(nsd+1,nenmax,ngaussmax),s(nee,nee)
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
      integer l
      double precision shd(4,nenmax,ngaussmax),b(6,3*nenmax),dtmp(6,6)
      double precision shbar(4,nenmax),db(6,3*nenmax),det(ngaussmax)
      double precision vol
c
c...form shape functions for each integration point
c
cdebug      write(6,*) "Hello from stiff_ss_f!"
c
      call getshape(xl,sh,shj,shd,shbar,det,gauss,val,iel,nen,nsd,
     & ngauss,ierr)
      if(ierr.ne.0) return
c
c...construct b matrix, then form intermediate db=dmat*b, and finally
c   the stiffness (b)t*dmat*b multiplied by appropriate weight for
c   integral over element
c
      do l=1,ngauss
        call bmatrix(b,shd(1,1,l),shbar,nsd,ndof,nen,nstr)
        call getmat(dtmp,dmat(1,l),iddmat,nstr,nddmat)
        call dsymm("l","l",nstr,nee,det(l),dtmp,nstr,b,nstr,zero,db,
     &   nstr)
        call dgemm("t","n",nee,nee,nstr,one,b,nstr,db,nstr,one,s,nee)
      end do
      return
      end
c
c version
c $Id: stiff_ss.f,v 1.1 2004/06/15 17:22:42 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
