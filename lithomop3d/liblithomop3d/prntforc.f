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
      subroutine prntforc(n,p,ien,nen,ndof,idout,kw)
c
c...prints local force vector for debugging
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer n,nen,ndof,idout,kw
      integer ien(nen)
      double precision p(ndof*nen)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  defined constants
c
      character*4 head(8)
      data head/8*'node'/
c
c...  local variables
c
      integer i,j
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
cdebug      write(6,*) "Hello from prntforc_f!"
c
      if(idout.lt.2) return
      write(kw,1000) n,(head(i),i,i=1,nen)
      write(kw,2000) (ien(i),i=1,nen)
      do j=1,ndof
        write(kw,3000) labeld(j),(p(j+(i-1)*ndof),i=1,nen)
      end do
1000  format(/,'local forces in element# ',i7/
     & 6x,'local   ',8(a4,i2,6x))
2000  format(6x,'global#  ',8(i7,7x))
3000  format(1x,a4,4x,8(1pe12.3))
      return
      end
c
c version
c $Id: prntforc.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
