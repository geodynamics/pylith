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
      subroutine printf(tfault,dfault,deltp,nfault,numfn,idout,
     & idsk,kw,kp)
c
c...prints out split node displacements
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c... subroutine arguments
c
      integer numfn,idout,idsk,kw,kp
      integer nfault(3,numfn)
      double precision tfault(ndof,numfn),dfault(ndof,numfn),deltp
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer npage,n,i,j
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
cdebug      write(6,*) "Hello from printf_f!"
c
      if(numfn.eq.0) return
      npage=50
      do n=1,numfn
        if((n.eq.1.or.mod(n,npage).eq.0).and.idout.gt.1) then
          write(kw,1000) (labeld(i),i=1,ndof)
          write(kw,*) ' '
        end if
        if(idout.gt.1) write(kw,2000) nfault(1,n),nfault(2,n),
     &   (tfault(i,n),i=1,ndof)
        if(deltp.eq.zero) then
          if(idsk.eq.0) write(kp,3000) (tfault(i,n),i=1,ndof)
        else
          if(idsk.eq.0) write(kp,3000) (tfault(i,n),i=1,ndof),
     &     (dfault(j,n)/deltp,j=1,ndof)
        end if
      end do
      if(idsk.eq.1) then
        write(kp) tfault
        if(deltp.ne.zero) then
          write(kp) deltp
          write(kp) dfault
        end if
      end if
1000  format(///1x,'s p l i t   n o d e   d i s p l a c e m e n t s '//
     &    1x,'element',4x,'node',15x,6(a4,14x))
2000  format(1x,i7,4x,i7,6x,6(1pe18.5))
3000  format(12e15.7)
      return
      end
c
c version
c $Id: printf.f,v 1.2 2004/07/07 20:27:07 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
