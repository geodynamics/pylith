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
      subroutine localf(ien,lmf,infiel,nconsz,numelt,infetype,nfault,
     & numfn)
c
c.......subroutine to localize nfault array for efficient computation
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nconsz,numelt,numfn
      integer ien(nconsz),lmf(nconsz),infiel(6,numelt)
      integer infetype(4,netypes),nfault(3,numfn)
c
c...  local variables
c
      integer i,j,iel,indien,ietype,nen,node
c
      call ifill(lmf,izero,nconsz)
      if(numfn.eq.izero) return
c
      do i=1,numfn
        iel=nfault(1,i)
        indien=infiel(1,iel)
        ietype=infiel(3,iel)
        nen=infetype(2,ietype)
        node=nfault(2,i)
        do j=indien,indien+nen-1
          if(node.eq.ien(j)) lmf(j)=i
        end do
      end do
      return
      end
c
c version
c $Id: localf.f,v 1.2 2004/07/07 18:27:14 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
