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
      subroutine getshapb(xl,sh,shj,shd,shbar,det,gauss,vol,iel,nen,nsd,
     & ngauss,ierr)
c
c...  Subroutine to compute shape functions and derivatives at
c     Gauss points.
c
c     This is a generic routine for any element type, and it uses
c     Bbar selective integration.
c
      include "implicit.inc"
c
c...  dimension parameters
c
      include "nshape.inc"
c
c...  subroutine arguments
c
      integer iel,nen,nsd,ngauss,ierr
      double precision xl(nsd,nen),sh(nsd+1,nenmax,ngaussmax)
      double precision shj(nsd+1,nenmax,ngaussmax)
      double precision shd(nsd+1,nenmax,ngaussmax),shbar(nsd+1,nenmax)
      double precision det(ngauss),gauss(nsd+1,ngaussmax),vol
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  local variables
c
cdebug      integer idb,jdb
      integer l
      double precision xs(3,3)
c
c...compute shape function derivatives over number of strain integration
c   points, and then compute average dilatational component.
c
cdebug      write(6,*) "Hello from getshapb_f!"
c
      vol=zero
      do l=1,ngauss
        call getjac(xl,xs,det(l),shj(1,1,l),nen,nsd,iel,ierr)
        if(ierr.ne.0) return
        call getder(det(l),sh(1,1,l),shd(1,1,l),xs,nen,nsd)
        det(l)=gauss(4,l)*det(l)
        vol=vol+det(l)
      end do
      call meansh(shbar,shd,vol,det,nsd,nen,ngauss)
      return
      end
c
c version
c $Id: getshapb.f,v 1.3 2004/06/17 18:09:43 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
