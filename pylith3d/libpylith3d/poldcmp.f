c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
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
      subroutine poldcmp(ee,wr,r,u)
c
c...routine to perform a polar decomposition on the deformation
c   gradient, given the linear strain and rotation matrices
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "rconsts.inc"
c
      integer idim
      parameter(idim=3)
c
c...  subroutine arguments
c
      double precision ee(6),wr(3),r(3,3),u(3,3)
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  local variables
c
      integer nrot
      double precision e1,e2,e3
      double precision c(3,3),x(3,3)
c
c...form inverse deformation gradient and inverse Cauchy deformation
c   tensor
c
cdebug      write(6,*) "Hello from poldcmp_f!"
c
      x(1,1)=one-ee(1)
      x(2,2)=one-ee(2)
      x(3,3)=one-ee(3)
      x(1,2)=-ee(4)+wr(1)
      x(2,1)=-ee(4)-wr(1)
      x(2,3)=-ee(5)+wr(3)
      x(3,2)=-ee(5)-wr(3)
      x(1,3)=-ee(6)+wr(2)
      x(3,1)=-ee(6)-wr(2)
      call dsyrk("u","n",idim,idim,one,x,idim,zero,c,idim)
      call symmet(c,idim)
c
c...perform eigenvalue decomposition
c
c******  need to use a LAPACK eigenvalue routine -- then symmetrization
c******  probably won't be necessary.
      call jacobi(c,idim,idim,ee,r,nrot)
c
c...compute stretch tensor and rotation matrix
c
      e1=zero
      e2=zero
      e3=zero
      if(ee(1).ne.zero) e1=one/sqrt(ee(1))
      if(ee(2).ne.zero) e2=one/sqrt(ee(2))
      if(ee(3).ne.zero) e3=one/sqrt(ee(3))
      u(1,1)=r(1,1)*r(1,1)*e1+r(1,2)*r(1,2)*e2+r(1,3)*r(1,3)*e3
      u(2,2)=r(2,1)*r(2,1)*e1+r(2,2)*r(2,2)*e2+r(2,3)*r(2,3)*e3
      u(3,3)=r(3,1)*r(3,1)*e1+r(3,2)*r(3,2)*e2+r(3,3)*r(3,3)*e3
      u(1,2)=r(1,1)*r(2,1)*e1+r(1,2)*r(2,2)*e2+r(1,3)*r(2,3)*e3
      u(1,3)=r(1,1)*r(3,1)*e1+r(1,2)*r(3,2)*e2+r(1,3)*r(3,3)*e3
      u(2,3)=r(2,1)*r(3,1)*e1+r(2,2)*r(3,2)*e2+r(2,3)*r(3,3)*e3
      u(2,1)=u(1,2)
      u(3,1)=u(1,3)
      u(3,2)=u(2,3)
      call dgemm("t","t",idim,idim,idim,one,x,idim,u,idim,zero,r,idim)
      return
      end
c
c version
c $Id: poldcmp.f,v 1.4 2005/04/08 00:41:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
