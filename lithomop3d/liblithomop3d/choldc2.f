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
      subroutine choldc2(a,n,np,p,ierr,errstrng)
c
c...  Subroutine to perform Cholesky decomposition.
c     Adapted from Numerical Recipes.
c...  Note:  See about replacing this routine with lapack routine
c     SPOFA (and then corresponding call to SPOSL).
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
      integer n,np,ierr
      character errstrng*(*)
      double precision a(np,np),p(n)
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  local variables
c
      integer i,j,k
      double precision sum
c
c*      write(6,*) "Hello from choldc2_f!"
c
      ierr=izero
c
      do i=1,n
        do j=i,n
          sum=a(i,j)
          do k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
          end do
          if(i.eq.j)then
            if(sum.le.zero)then
              ierr=110
              errstrng="choldc2"
              return
            end if
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          end if
        end do
      end do
      return
      end
c
c version
c $Id: choldc2.f,v 1.4 2004/06/21 20:08:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
