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
      subroutine bmatrixb(b,sh,shbar,nen)
c
c...computes the linear strain-displacement matrix b(nstr,ndof*nen)
c
c      This routine uses Hughes' B-bar formuation for dilatational
c      strains.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nen
      double precision b(nstr,ndof*nenmax),sh(nsd+1,nenmax)
      double precision shbar(nsd+1,nenmax)
c
c...  local variables
c
      integer k,i
      double precision dil1,dil2,dil3
c
cdebug      write(6,*) "Hello from bmatrixb_f!"
c
      k=1
      do i=1,nen
        dil1=third*(shbar(1,i)-sh(1,i))
        dil2=third*(shbar(2,i)-sh(2,i))
        dil3=third*(shbar(3,i)-sh(3,i))
        b(1,k  )=sh(1,i)+dil1
        b(1,k+1)=dil2
        b(1,k+2)=dil3
        b(2,k  )=dil1
        b(2,k+1)=sh(2,i)+dil2
        b(2,k+2)=dil3
        b(3,k  )=dil1
        b(3,k+1)=dil2
        b(3,k+2)=sh(3,i)+dil3
        b(4,k  )=sh(2,i)
        b(4,k+1)=sh(1,i)
        b(4,k+2)=zero
        b(5,k  )=zero
        b(5,k+1)=sh(3,i)
        b(5,k+2)=sh(2,i)
        b(6,k  )=sh(3,i)
        b(6,k+1)=zero
        b(6,k+2)=sh(1,i)
        k=k+3
      end do
      return
      end
c
c version
c $Id: bmatrixb.f,v 1.3 2004/08/12 01:06:26 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
