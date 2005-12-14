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
      subroutine plinhex(sh,gauss,nen,ngauss,intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a trilinear hexahedron.
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
      double precision r(8),s(8),t(8)
      data r/-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0/
      data s/-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0/
      data t/ 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0/
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i,l
      double precision rr,ss,tt
c
cdebug      integer idb
c
cdebug      write(6,*) "Hello from plinhex_f!"
c
c...  definitions
c
c
c...  Linear hex definition
c
      gauss(1,1)=zero
      gauss(2,1)=zero
      gauss(3,1)=zero
      gauss(4,1)=eight
      if(intord.ne.2) then
        do l=1,ngauss
          gauss(1,l)=r(l)*root3i
          gauss(2,l)=s(l)*root3i
          gauss(3,l)=t(l)*root3i
          gauss(4,l)=one
        end do
      end if
c
      do l=1,ngauss
        do i=1,nen
          rr=one+r(i)*gauss(1,l)
          ss=one+s(i)*gauss(2,l)
          tt=one+t(i)*gauss(3,l)
          sh(4,i,l)=eighth*rr*ss*tt
          sh(1,i,l)=eighth*r(i)*ss*tt
          sh(2,i,l)=eighth*s(i)*rr*tt
          sh(3,i,l)=eighth*t(i)*rr*ss
        end do
      end do
c
      return
      end
c
c version
c $Id: plinhex.f,v 1.4 2005/03/22 04:45:54 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
