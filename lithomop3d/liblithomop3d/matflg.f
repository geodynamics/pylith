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
      subroutine matflg(prop,nprop,ivisc,iplas,imhist)
c
c...subroutine to set viscous and plastic solution flags, depending on
c   the values of the material properties
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nprop,ivisc,iplas,imhist
      double precision prop(nprop)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
      double precision big
      parameter(big=1.0d30)
c
c...  local variables
c
      integer i
c
c*      write(6,*) "Hello from matflg_f!"
c
      if(prop(4).lt.big) ivisc=ione
      if(prop(7).lt.big) iplas=ione
      do i=1,nprop
        if(prop(i).lt.zero) imhist=ione
      end do
      return
      end
c
c version
c $Id: matflg.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
