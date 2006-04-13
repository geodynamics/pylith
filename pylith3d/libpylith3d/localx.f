c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
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
      subroutine localx(idx,numnp,ien,lmx,numelv,nslip,numslp,nen)
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
      integer numnp,numelv,numslp,nen
      integer idx(ndof,numnp),ien(nen,numelv),lmx(ndof,nen,numelv)
      integer nslip(nsdim,numslp)
c
c...  intrinsic functions
c
      intrinsic sign
c
c...  local variables
c
      integer i,j,k,ielg,node
c
cdebug      integer idb
c
cdebug      write(6,*) "Hello from localx_f!"
c
      call ifill(lmx,izero,ndof*nen*numelv)
      if(numslp.eq.izero) return
c
      do i=1,numslp
        ielg=nslip(1,i)
        node=nslip(2,i)
        do j=1,nen
          if(node.eq.ien(j,ielg)) then
            do k=1,ndof
              lmx(k,j,ielg)=sign(idx(k,node),nslip(2+k,i))
            end do
          end if
cdebug          write(6,*) "i,j,node,ielg,lmx:",
cdebug     &     i,j,node,ielg,(lmx(idb,j,ielg),idb=1,ndof)
        end do
      end do
      return
      end
c
c version
c $Id: localx.f,v 1.6 2005/04/08 00:31:37 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
