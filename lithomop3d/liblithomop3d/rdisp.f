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
      subroutine rdisp(d,skew,numnp)
c
c...rotates displacement vector to global coordinate directions
c***  Note:  this could probably be done much more efficiently
c     using a skew array of dimension numrot in conjunction with
c     an index array (of the same dimension).
c     Also note:  an additional index array could keep track of which
c     skew BC correspond to slippery nodes (for the case of automatic
c     skew computation).
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
      integer numnp
      double precision d(ndof,numnp),skew(nskdim,numnp)
c
c...  local variables
c
      integer i
      double precision rot(3,3),dtemp(3)
c
      do i=1,numnp
        if((skew(1,i).ne.zero).and.(skew(nskdim,i).ne.zero)) then
          call formrt(skew(1,i),rot,nskdim)
	  call dcopy(ndof,d(1,i),ione,dtemp,ione)
	  call dgemv("n",ndof,ndof,one,rot,ithree,dtemp,ione,zero,
     &     d(1,i),ione)
        end if
      end do
      return
      end
c
c version
c $Id: rdisp.f,v 1.2 2004/07/01 21:03:17 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
