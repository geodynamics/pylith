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
      subroutine stiffldql(sh,stn,s,det,ngauss,ndof,nstr,nen,nee)
c
c...subroutine to compute additional stiffness contribution for
c   large deformations and add it to regular stiffness
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ngauss,ndof,nstr,nen,nee
      double precision sh(4,8,ngauss),stn(nstr,ngauss),s(nee,nee)
      double precision det(ngauss)
c
c...  defined constants
c
      include "rconsts.inc"
c
      integer idim
      parameter(idim=9)
c
c...  local variables
c
      integer l,i
      double precision bn(idim,24),sb(idim,24),stnm(idim,idim)
c
c...construct bn matrix, then form intermediate sb=stn*bn, and finally
c   the stiffness (bn)t*stn*bn multiplied by appropriate weight for
c   integral over element
c
c*      write(6,*) "Hello from stiffldql_f!"
c
      call fill(stnm,zero,idim*idim)
      do l=1,ngauss
        call bnmatrxql(bn,sh(1,1,l))
        do i=1,3
          stnm(3*(i-1)+1,3*(i-1)+1)=stn(1,l)
          stnm(3*(i-1)+2,3*(i-1)+2)=stn(2,l)
          stnm(3*(i-1)+3,3*(i-1)+3)=stn(3,l)
          stnm(3*(i-1)+1,3*(i-1)+2)=stn(4,l)
          stnm(3*(i-1)+2,3*(i-1)+1)=stn(4,l)
          stnm(3*(i-1)+2,3*(i-1)+3)=stn(5,l)
          stnm(3*(i-1)+3,3*(i-1)+2)=stn(5,l)
          stnm(3*(i-1)+1,3*(i-1)+3)=stn(6,l)
          stnm(3*(i-1)+3,3*(i-1)+1)=stn(6,l)
        end do
        call dsymm("l","l",idim,nee,det(l),stnm,idim,bn,idim,zero,sb,
     &   idim)
        call dgemm("t","n",nee,nee,idim,one,bn,idim,sb,idim,zero,s,nee)
      end do
      return
      end
c
c version
c $Id: stiffldql.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
