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
      subroutine getder(det,sh,shd,xs,nen,nsd)
c
c...  subroutine to compute shape function derivatives in
c     global coordinates
c
c        sh(1,nen),sh(2,nen),sh(3,nen)    = r, s, and t derivatives
c                                           of shape functions
c        shd(1,nen),shd(2,nen),shd(3,nen) = x, y, and z derivatives
c                                           of shape functions
c        xs(nsd,nsd)                      = jacobian matrix
c        det                              = jacobian matrix determinant
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nen,nsd
      double precision det,sh(nsd+1,nen),shd(nsd+1,nen),xs(nsd,nsd)
c
c...  defined constants
c
      include "rconsts.inc"
c
c...  local variables
c
      integer i
      double precision a11,a12,a13,a21,a22,a23,a31,a32,a33,detinv
c
cdebug      write(6,*) "Hello from getder_f!"
c
      detinv=one/det
c
c...transform natural derivatives to (x,y,z) derivatives using an
c   explicit 3x3 matrix inversion routine
c
      a11 = (xs(2,2)*xs(3,3))-(xs(2,3)*xs(3,2))
      a12 =-(xs(1,2)*xs(3,3))+(xs(1,3)*xs(3,2))
      a13 = (xs(1,2)*xs(2,3))-(xs(1,3)*xs(2,2))
      a21 =-(xs(2,1)*xs(3,3))+(xs(2,3)*xs(3,1))
      a22 = (xs(1,1)*xs(3,3))-(xs(1,3)*xs(3,1))
      a23 =-(xs(1,1)*xs(2,3))+(xs(1,3)*xs(2,1))
      a31 = (xs(2,1)*xs(3,2))-(xs(2,2)*xs(3,1))
      a32 =-(xs(1,1)*xs(3,2))+(xs(1,2)*xs(3,1))
      a33 = (xs(1,1)*xs(2,2))-(xs(1,2)*xs(2,1))
      detinv=one/det
      do i=1,nen
        shd(1,i)=detinv*(a11*sh(1,i)+a12*sh(2,i)+a13*sh(3,i))
        shd(2,i)=detinv*(a21*sh(1,i)+a22*sh(2,i)+a23*sh(3,i))
        shd(3,i)=detinv*(a31*sh(1,i)+a32*sh(2,i)+a33*sh(3,i))
        shd(4,i)=sh(4,i)
      end do
c
      return
      end
c
c version
c $Id: getder.f,v 1.2 2004/06/15 19:44:47 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
