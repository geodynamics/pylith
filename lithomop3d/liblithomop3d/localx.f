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
      subroutine localx(idx,ien,lmx,nslip,nen,ndof,numslp,numel,numnp,
     & nsdim)
c
c......subroutine to localize idx array and transfer sign information
c       to lmx array
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nen,ndof,numslp,numel,numnp,nsdim
      integer idx(ndof,numnp),ien(nen,numel),lmx(ndof,nen,numel)
      integer nslip(nsdim,numslp)
c
c...  defined constants
c
      include "nconsts.inc"
c
c...  intrinsic functions
c
      intrinsic sign
c
c...  local variables
c
      integer i,j,k,node,nelem
c
      call ifill(lmx,izero,ndof*nen*numel)
      if(numslp.eq.0) return
c
      do i=1,numslp
        node=nslip(2,i)
        nelem=nslip(1,i)
        do j=1,nen
          if(node.eq.ien(j,nelem)) then
            do k=1,ndof
              lmx(k,j,nelem)=sign(idx(k,node),nslip(2+k,i))
            end do
          end if
        end do
      end do
      return
      end
c
c version
c $Id: localx.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
