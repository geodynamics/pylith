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
      subroutine gspre(alnz,bres,z,rtz,ja,neq,nnz)

c...  subroutine to perform symmetrized Gauss-Seidel preconditioning
c     This routine returns z=B-inverse*bres, where B is an approximation
c     to the matrix alnz.  The routine also returns the dot product of
c     z and bres.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer neq,nnz
      integer ja(nnz)
      double precision alnz(nnz),bres(neq),z(neq),rtz
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  local variables
c
      integer i,j
      double precision sum,dsr
c
c...  perform forward reduction
c
      call dcopy(neq,bres,ione,z,ione)
      do i=1,neq
        sum=z(i)
        do j=ja(i),ja(i+1)-1
          if(ja(j).gt.i) go to 10
          sum=sum-z(ja(j))*alnz(j)/sqrt(alnz(ja(j)))
        end do
10      continue
        z(i)=sum/sqrt(alnz(i))
      end do
c
c...  perform back substitution
c
      rtz=zero
      do i=neq,1,-1
        dsr=one/sqrt(alnz(i))
        sum=z(i)
        do j=ja(i+1)-1,ja(i),-1
          if(ja(j).lt.i) go to 20
          sum=sum-z(ja(j))*alnz(j)*dsr
        end do
20      continue
        z(i)=sum*dsr
        rtz=rtz+bres(i)*z(i)
      end do
      return
      end
c
c version
c $Id: gspre.f,v 1.2 2004/07/07 14:17:41 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
