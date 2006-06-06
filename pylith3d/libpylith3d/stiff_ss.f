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
     & xl,iddmat,                                                       ! global
     & dmat,ien,iel,                                                    ! elemnt
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & s,                                                               ! stiff
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...computes the local stiffness matrix at the given integration points.
c   k=(b)t*d*b.
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
      integer iel,nen,ngauss,nee,ierr
      integer iddmat(nstr,nstr),ien(nen)
      character errstrng*(*)
      double precision xl(nsd,nen),dmat(nddmat,ngauss)
      double precision gauss(nsd+1,ngauss),sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss),s(nee,nee)
c
c...  external routines
c
      external getshape,bmatrix
c
c...  local variables
c
      integer l
      double precision shd(4,nenmax,ngaussmax),b(6,3*nenmax)
      double precision shbar(4,nenmax),db(6,3*nenmax),det(ngaussmax)
      double precision dtmp(6,6)
      double precision vol
c
cdebug      integer idb,jdb
c
c...form shape functions for each integration point
c
cdebug      write(6,*) "Hello from stiff_ss_f!"
c
      call getshape(xl,sh,shj,shd,shbar,det,gauss,vol,iel,nen,ngauss,
     & ierr,errstrng)
      if(ierr.ne.izero) return
c
c...construct b matrix, then form intermediate db=dmat*b, and finally
c   the stiffness (b)t*dmat*b multiplied by appropriate weight for
c   integral over element
c
      do l=1,ngauss
        call bmatrix(b,shd(1,1,l),shbar,nen)
        call getmat(dtmp,dmat(1,l),iddmat,nstr,nddmat)
        call dsymm("l","l",nstr,nee,det(l),dtmp,nstr,b,nstr,zero,db,
     &    nstr)
c*        call dspmv("u",nstr,det(l),dmat(1,l),b,ione,zero,db,ione)
        call dgemm("t","n",nee,nee,nstr,one,b,nstr,db,nstr,one,s,nee)
cdebug        if(iel.eq.8) then
cdebug          write(6,*) "det:",det(l)
cdebug          do idb=1,nee
cdebug            write(6,*) "b:",(b(jdb,idb),jdb=1,nstr)
cdebug          end do
cdebug          write(6,*) "dmat:",(dmat(jdb,l),jdb=1,nddmat)
cdebug        end if
      end do
cdebug      do idb=1,nee
cdebug        write(6,*) "idb,s(idb,idb):",idb,s(idb,idb)
cdebug      end do
      return
      end
c
c version
c $Id: stiff_ss.f,v 1.6 2005/04/01 23:31:35 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
