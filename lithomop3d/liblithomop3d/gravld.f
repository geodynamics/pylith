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
      subroutine gravld(p,grav,xl,iel,nen,dens,gauss,shj,ngauss,ierr,
     & errstrng)
c
c...computes the local forces due to gravitational acceleration
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
      integer iel,nen,ngauss,ierr
      character errstrng*(*)
      double precision p(ndof*nen),grav(ndof),xl(nsd,nen),dens
      double precision gauss(nsd+1,ngauss)
      double precision shj(nsd+1,nen,ngauss)
c
c...  local variables
c
      integer l,j
      double precision xs(3,3),det,rl1,rl2,rl3
c
cdebug      write(6,*) "Hello from gravld_f!"
c
      if((grav(1)*grav(1)+grav(2)*grav(2)+grav(3)*grav(3)).eq.zero)
     & return
      if(dens.eq.zero) return
      do l=1,ngauss
        call getjac(xl,xs,det,shj(1,1,l),nen,iel,ierr,errstrng)
        if(ierr.ne.izero) return
        det=det*gauss(4,l)*dens
        rl1=det*grav(1)
        rl2=det*grav(2)
        rl3=det*grav(3)
        do j=1,nen
          p(3*j-2)=p(3*j-2)+rl1*shj(4,j,l)
          p(3*j-1)=p(3*j-1)+rl2*shj(4,j,l)
          p(3*j  )=p(3*j  )+rl3*shj(4,j,l)
        end do
      end do
      return
      end
c
c version
c $Id: gravld.f,v 1.5 2005/03/21 22:13:15 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
