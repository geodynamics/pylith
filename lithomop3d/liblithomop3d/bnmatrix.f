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
      subroutine bnmatrix(bn,sh,nen)
c
c...computes the nonlinear strain-displacement matrix
c   bn(ndof*ndof,ndof*nen) used in computing the initial stress
c   stiffness for 3 dimensional problems.
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
      double precision bn(ndof*ndof,ndof*nenmax),sh(nsd+1,nenmax)
c
c...  local variables
c
      integer k,i
c
c*      write(6,*) "Hello from bnmatrix_f!"
c
      k=1
      do i=1,nen
        bn(1,k  )=sh(1,i)
        bn(1,k+1)=zero
        bn(1,k+2)=zero
        bn(2,k  )=sh(2,i)
        bn(2,k+1)=zero
        bn(2,k+2)=zero
        bn(3,k  )=sh(3,i)
        bn(3,k+1)=zero
        bn(3,k+2)=zero
        bn(4,k  )=zero
        bn(4,k+1)=sh(1,i)
        bn(4,k+2)=zero
        bn(5,k  )=zero
        bn(5,k+1)=sh(2,i)
        bn(5,k+2)=zero
        bn(6,k  )=zero
        bn(6,k+1)=sh(3,i)
        bn(6,k+2)=zero
        bn(7,k  )=zero
        bn(7,k+1)=zero
        bn(7,k+2)=sh(1,i)
        bn(8,k  )=zero
        bn(8,k+1)=zero
        bn(8,k+2)=sh(2,i)
        bn(9,k  )=zero
        bn(9,k+1)=zero
        bn(9,k+2)=sh(3,i)
        k=k+3
      end do
      return
      end
c
c version
c $Id: bnmatrix.f,v 1.1 2004/06/18 15:14:02 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
