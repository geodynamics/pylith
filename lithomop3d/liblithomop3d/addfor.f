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
      subroutine addfor(b,p,lm,lmx,neq,nee)
c
c...program to add effective force, including slippery node degrees
c   of freedom
c
      include "implicit.inc"
c
c... parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer neq,nee
      integer lm(nee),lmx(nee)
      double precision b(neq),p(nee)
c
c...  intrinsic functions
c
      intrinsic abs,sign,dble
c
c...  local variables
c
      integer j,k,l
      double precision sgn
c
cdebug      write(6,*) "Hello from addfor_f!"
c
      do j=1,nee
        k=lm(j)
        l=abs(lmx(j))
        sgn=sign(one,dble(lmx(j)))
cdebug        write(6,*) j,k,l,sgn
cdebug        write(6,*) j,p(j)
        if(k.ne.izero) b(k)=b(k)+p(j)
        if(l.ne.izero) b(l)=b(l)+p(j)*sgn
cdebug        if(k.ne.izero) write(6,*) k,b(k)
      end do
      return
      end
c
c version
c $Id: addfor.f,v 1.4 2004/08/02 21:02:59 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
