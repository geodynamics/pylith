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
      subroutine infcmp(ietypei,ietype,inf)
c
c...  routine that converts primitive element type plus infinite element
c     info to global element type.  The assumption at present is that
c     things are set up the same way as in the preshape routine.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer ietypei,ietype,inf
c
c...  local constants
c
      integer idiv(3)
      data idiv/1,10,100/
c
c...  local variables
c
      integer io(3),ind
c
cdebug      write(6,*) "Hello from infcmp_f!"
c
c
c...  simplest case:  no infinite elements
c
      if(inf.eq.izero) then
        if(ietypei.eq.ione) ietype=ione
        if(ietypei.gt.ione) ietype=ietypei+26
        if(ietypei.gt.isix) ietype=ietypei+52
      else if(ietypei.eq.ione) then
        io(3)=inf/idiv(3)
        io(2)=(inf-io(3)*idiv(3))/idiv(2)
        io(1)=inf-io(3)*idiv(3)-io(2)*idiv(2)
        ind=io(1)+ithree*io(2)+inine*io(3)
        ietype=ietypei+ind
      else if(ietypei.eq.isix) then
        io(3)=inf/idiv(3)
        io(2)=(inf-io(3)*idiv(3))/idiv(2)
        io(1)=inf-io(3)*idiv(3)-io(2)*idiv(2)
        ind=io(1)+ithree*io(2)+inine*io(3)
        ietype=ietypei+ind+26
      end if
      return
      end
c
c version
c $Id: infcmp.f,v 1.2 2004/08/12 01:30:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
