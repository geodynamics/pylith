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
      subroutine rpforc(p,skew,ien,numnp,nen)
c
c...rotates local force vector into skewed coordinate system
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numnp,nen
      integer ien(nen)
      double precision p(ndof*nen),skew(nskdim,numnp)
c
c...  local variables
c
      integer i,k,ll
      double precision rot(3,3),ptemp(3)
c
      do i=1,nen
        k=ien(i)
        if((skew(1,k).ne.zero).and.(skew(itwo,k).ne.zero)) then
          ll=ndof*(i-1)
          call formrt(skew(1,k),rot)
	  call dcopy(ndof,p(ll+1),ione,ptemp,ione)
	  call dgemv("t",ndof,ndof,one,rot,ithree,ptemp,ione,zero,
     &     p(ll+1),ione)
        end if
      end do
      return
      end
c
c version
c $Id: rpforc.f,v 1.2 2004/07/01 19:59:55 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
