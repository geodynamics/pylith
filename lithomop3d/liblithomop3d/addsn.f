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
      subroutine addsn(dl,dx,ien,lmx,nen,numnp)
c
c...adds displacements due to slip across internal interfaces to
c   the local displacement vectors
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
      integer nen,numnp
      integer ien(nen),lmx(ndof,nen)
      double precision dl(ndof,nen),dx(ndof,numnp)
c
c...  intrinsic functions
c
      intrinsic sign,dble
c
c...  local variables
c
      integer i,j,k
      double precision sgn
      logical*4 slip
c
      do i=1,nen
        k=ien(i)
        slip=.false.
        do j=1,ndof
          if(lmx(j,i).ne.izero) then
            slip=.true.
            sgn=sign(one,dble(lmx(j,i)))
          end if
        end do
        if(slip) then
          do j=1,ndof
            dl(j,i)=dl(j,i)+sgn*dx(j,k)
          end do
        end if
      end do
      return
      end
c
c version
c $Id: addsn.f,v 1.3 2004/07/06 20:12:00 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
