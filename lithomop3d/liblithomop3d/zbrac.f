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
      subroutine zbrac(x1,x2,succes,ae,bs,c,ds,dl1,dl2,t1,t2,emhu,
     & anpwr,efstsi,gam,dlam,deltp,alfap,iopd,plas)
c
c...subroutine to bracket the root of the effective stress function
c     adapted from Numerical Recipes
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iopd
      double precision x1,x2,ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,efstsi
      double precision gam,dlam,deltp,alfap
      logical succes,plas
c
c...  defined constants
c
      include "rconsts.inc"
c
      integer ntry
      double precision factor
      parameter (factor=1.6d0,ntry=100)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer j
      double precision f1,df1,f2,df2
c
cdebug2      write(6,*) "Hello from zbrac_f!"
c
      if(x1.eq.x2)pause 'you have to guess an initial range'
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x1,efstsi,gam,dlam,
     & f1,df1,deltp,alfap,iopd,plas)
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x2,efstsi,gam,dlam,
     & f2,df2,deltp,alfap,iopd,plas)
cdebug2      write(6,*) "Point 1 in zbrac"
      succes=.true.
      do 11 j=1,ntry
        if(f1*f2.lt.zero)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+factor*(x1-x2)
          if(x1.lt.zero) x1=zero
          call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x1,efstsi,gam,
     &     dlam,f1,df1,deltp,alfap,iopd,plas)
        else
          x2=x2+factor*(x2-x1)
          call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x2,efstsi,gam,
     &     dlam,f2,df2,deltp,alfap,iopd,plas)
        end if
11    continue
      succes=.false.
      return
      end
c
c version
c $Id: zbrac.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
