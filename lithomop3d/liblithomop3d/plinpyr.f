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
      subroutine plinpyr(sh,shj,gauss,ngauss,nen,nsd,nenmax,ngaussmax,
     & intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a linear pyramid.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nsd,nenmax,ngaussmax,intord
      integer ngauss,nen
      double precision sh(nsd+1,nenmax,ngaussmax)
      double precision shj(nsd+1,nenmax,ngaussmax)
      double precision gauss(nsd+1,ngaussmax)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
      double precision c1,c2,c3,c4,c5,c6
      parameter(c1 = 128.0d0,
     &          c2 =  27.0d0,
     &          c3 =  15.0d0,
     &          c4 =  81.0d0,
     &          c5 = 100.0d0,
     &          c6 = 125.0d0)
c
      double precision r(5),s(5),t(5)
      data r/-1d0, 1d0, 1d0,-1d0, 0d0/
      data s/-1d0,-1d0, 1d0, 1d0, 0d0/
      data t/-1d0,-1d0,-1d0,-1d0, 1d0/
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i,l,nshsize,ngssize
      double precision rr,ss,tt,g1,w1
c
c...  definitions
c
      nshsize=(nsd+1)*nenmax*ngaussmax
      ngssize=(nsd+1)*ngaussmax
c
c...  Linear wedge definition
c
      nen=ifive
      ngauss=ione
      gauss(1,1)=zero
      gauss(2,1)=zero
      gauss(3,1)=-half
      gauss(4,1)=c1/c2
      if(intord.ne.2) then
        ngauss=ifive
        g1=eight*sqrt(two/c3)/five
        w1=c4/c5
        do l=1,ngauss
          gauss(1,l)=r(l)*g1
          gauss(2,l)=s(l)*g1
          gauss(3,l)=-two*third
          gauss(4,l)=w1
        end do
        gauss(3,5)=two/five
        gauss(4,5)=c6/c2
      end if
      do l=1,ngauss
        do i=1,nen
          rr=half*(one+r(i)*gauss(1,l))
          ss=half*(one+s(i)*gauss(2,l))
          tt=half*(one-t(i)*gauss(3,l))
          sh(4,i,l)=rr*ss*tt
          sh(1,i,l)=r(i)*ss*tt
          sh(2,i,l)=s(i)*rr*tt
          sh(3,i,l)=t(i)*rr*ss
        end do
      end do
      call dcopy(nshsize,sh,ione,shj,ione)
c
      return
      end
c
c version
c $Id: plinpyr.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
