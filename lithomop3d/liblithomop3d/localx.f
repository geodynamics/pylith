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
      subroutine localx(idx,numnp,ien,lmx,infiel,nconsz,numelt,infetype,
     & nslip,numslp)
c
c......subroutine to localize idx array and transfer sign information
c       to lmx array
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
      integer numnp,nconsz,numelt,numslp
      integer idx(ndof,numnp),ien(nconsz),lmx(ndof,nconsz)
      integer infiel(7,numelt),infetype(4,netypes),nslip(nsdim,numslp)
c
c...  intrinsic functions
c
      intrinsic sign
c
c...  local variables
c
      integer i,j,k,iel,indien,ietype,nen,node
c
      call ifill(lmx,izero,ndof*nconsz)
      if(numslp.eq.izero) return
c
      do i=1,numslp
        iel=nslip(1,i)
        indien=infiel(1,iel)
        ietype=infiel(3,iel)
        nen=infetype(2,ietype)
        node=nslip(2,i)
        do j=indien,indien+nen-1
          if(node.eq.ien(j)) then
            do k=1,ndof
              lmx(k,j)=sign(idx(k,node),nslip(2+k,i))
            end do
          end if
        end do
      end do
      return
      end
c
c version
c $Id: localx.f,v 1.3 2005/02/24 00:03:56 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
