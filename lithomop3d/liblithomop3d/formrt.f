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
      subroutine formrt(skew,rot,nskdim)
c
c...constructs a three dimensional rotation matrix rot from the euler
c   angles contained in the array skew.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nskdim
      double precision skew(nskdim),rot(3,3)
c
c...  defined constants
c
      include "rconsts.inc"
c
c...  intrinsic functions
c
      intrinsic sin,cos
c
c...  local variables
c
      double precision sina,cosa,sinb,cosb
c
c*      write(6,*) "Hello from formrt_f!"
c
      sina=sin(skew(1))
      cosa=cos(skew(1))
      sinb=sin(skew(2))
      cosb=cos(skew(2))
      rot(1,1)=cosa*cosb
      rot(2,1)=sina*cosb
      rot(3,1)=-sinb
      rot(1,2)=-sina
      rot(2,2)=cosa
      rot(3,2)=zero
      rot(1,3)=cosa*sinb
      rot(2,3)=sina*sinb
      rot(3,3)=cosb
      return
      end
c
c version
c $Id: formrt.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
