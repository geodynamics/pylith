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
      subroutine rsplit(nfault,dfault,skew,ndof,numfn,numnp,nskdim)
c
c...rotates split node displacements to global coordinate directions
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ndof,numfn,numnp,nskdim
      integer nfault(3,numfn)
      double precision dfault(ndof,numfn),skew(nskdim,numnp)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  local variables
c
      integer i,node
      double precision rot(9),ftemp(3)
c
      do i=1,numfn
        node=nfault(2,i)
        if((skew(1,node).ne.zero).and.(skew(nskdim,node).ne.zero)) then
          call formrt(skew(1,node),rot,nskdim)
          call dcopy(ndof,dfault(1,i),ione,ftemp,ione)
          call dgemv("n",ndof,ndof,one,rot,ithree,ftemp,ione,zero,
     &     dfault(1,i),ione)
        end if
      end do
      return
      end
c
c version
c $Id: rsplit.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
