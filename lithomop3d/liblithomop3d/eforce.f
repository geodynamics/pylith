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
      subroutine eforce(xl,sh,shj,det,gauss,evp,p,iel,nen,ngauss,
     & getshape,bmatrix,ierr,errstrng)
c
c...this subroutine computes the effective forces at each node
c   within an element: p=(b)t*evp, where p is the local force vector,
c   b is the strain-displacement matrix, and evp is the stress
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "rconsts.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer iel,nen,ierr
      character errstrng*(*)
      double precision xl(nsd,nen),sh(nsd+1,nenmax,ngaussmax)
      double precision shj(nsd+1,nenmax,ngaussmax),det(ngaussmax)
      double precision gauss(nsd+1,ngaussmax),evp(nstr,ngauss)
      double precision p(ndof*nen)
      external getshape,bmatrix
c
c...  local variables
c
cdebug      integer idb,jdb
      integer l,nee
      double precision shd(4,nenmax,ngaussmax),shbar(4,20),b(6,60),vol
c
cdebug      write(6,*) "Hello from eforce_f!"
c
      nee=ndof*nen
      call getshape(xl,sh,shj,shd,shbar,det,gauss,vol,iel,nen,ngauss,
     & ierr,errstrng)
      if(ierr.ne.izero) return
c
      do l=1,ngauss
        call bmatrix(b,sh(1,1,l),shbar,nen)
        call dgemv("t",nstr,nee,det(l),b,nstr,evp(1,l),ione,one,p,ione)
      end do
      return
      end
c
c version
c $Id: eforce.f,v 1.3 2004/06/21 20:11:19 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
