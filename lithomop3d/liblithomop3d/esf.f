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
      subroutine esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,
     & efsts,efstsi,gam,dlam,val,dval,deltp,alfap,iopd,plas)
c
c...program to compute the effective stress function and its
c   derivative for a given effective strain and initial stress state
c
c**** Note that this routine is currently set up for the previous
c     method of determining material behavior (only one material
c     type with the option of including viscous or plastic behavior).
c     In the next version of the code, there should be a different
c     routine like this for every material type (except for purely
c     elastic, which won't need it).  That way, it won't be necessary
c     to perform switches based on whether viscous or plastic behavior
c     is being used.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iopd
      double precision ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,efsts,efstsi
      double precision gam,dlam,val,dval,deltp,alfap
      logical plas
c
c...  defined constants
c
      include "rconsts.inc"
c
c...  local variables
c
      double precision strtau,a,dgam,ddlam,ddl1,ddl2
c
c*      write(6,*) "Hello from esf_f!"
      dlam=dl1*efsts+dl2
      if(dlam.le.zero.or.(.not.plas)) dlam=zero
c*      if(dlam.le.zero) dlam=zero
      strtau=(one-alfap)*efstsi+alfap*efsts
      gam=half*(strtau/emhu)**(anpwr-one)/emhu
      a=ae+t2*gam
      if(efsts.ne.zero) a=a+dlam/(t1*efsts)
      val=a*a*efsts*efsts-bs*gam*gam+c*gam-ds
      if(iopd.eq.1) return
      dgam=alfap*(anpwr-one)*gam/strtau
      ddlam=dl1
      if(dlam.eq.zero) ddlam=zero
      ddl1=zero
      ddl2=zero
      if(efsts.ne.zero) then
        ddl1=ddlam/(t1*efsts)
        ddl2=dlam/(t1*efsts*efsts)
      end if
      dval=two*a*a*efsts-two*bs*gam*dgam+c*dgam+
     & two*a*efsts*efsts*(ddl1-ddl2+deltp*dgam*alfap)
      return
      end
c
c version
c $Id: esf.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
