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
      subroutine cholsl(a,n,np,p,b,x)
c
c...  routine to solve a square inverse problem given a
c     Cholesky-factored matrix a.
c     Taken from Numerical Recipes
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer n,np
      double precision a(np,np),b(n),p(n),x(n)
c
c...  local variables
c
      integer i,k
      double precision sum
c
cdebug      write(6,*) "Hello from cholsl_f!"
c
      do i=1,n
        sum=b(i)
        do k=i-1,1,-1
          sum=sum-a(i,k)*x(k)
        end do
        x(i)=sum/p(i)
      end do
      do i=n,1,-1
        sum=x(i)
        do k=i+1,n
          sum=sum-a(k,i)*x(k)
        end do
        x(i)=sum/p(i)
      end do
      return
      END
c
c version
c $Id: cholsl.f,v 1.2 2004/08/12 01:12:46 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
