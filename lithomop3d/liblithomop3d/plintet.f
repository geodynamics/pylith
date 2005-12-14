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
      subroutine plintet(sh,gauss,nen,ngauss,intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a linear tetrahedron.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nen,ngauss,intord
      double precision sh(nsd+1,nen,ngauss)
      double precision gauss(nsd+1,ngauss)
c
c...  local constants
c
      double precision r(4),s(4),t(4),u(4)
      data u/ 1d0, 0d0, 0d0, 0d0/
      data r/ 0d0, 1d0, 0d0, 0d0/
      data s/ 0d0, 0d0, 1d0, 0d0/
      data t/ 0d0, 0d0, 0d0, 1d0/
c
c...  intrinsic functions
c
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i
      double precision rr,ss,tt,uu
      double precision tetvol
c
c...  definitions
c
      tetvol=sixth
c
c...  Linear hex definition
c     One-point integration is used in all cases.
c
      gauss(1,1)=fourth
      gauss(2,1)=fourth
      gauss(3,1)=fourth
      gauss(4,1)=tetvol
c
      do i=1,nen
        rr=r(i)*gauss(1,1)
        ss=s(i)*gauss(2,1)
        tt=t(i)*gauss(3,1)
        uu=u(i)*gauss(3,1)
        sh(4,i,1)=rr+ss+tt+uu
        sh(1,i,1)=r(i)-u(i)
        sh(2,i,1)=s(i)-u(i)
        sh(3,i,1)=t(i)-u(i)
      end do
c
      return
      end
c
c version
c $Id: plintet.f,v 1.4 2005/03/22 04:45:55 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
