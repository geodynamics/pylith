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
      subroutine plinwedge(sh,shj,gauss,infetype,intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a linear wedge.
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
      integer intord
      integer infetype(4)
      double precision sh(nsd+1,nenmax,ngaussmax)
      double precision shj(nsd+1,nenmax,ngaussmax)
      double precision gauss(nsd+1,ngaussmax)
c
c...  local constants
c
      double precision r(6),s(6),t(6),u(6)
      data r/ 1d0, 0d0, 0d0, 1d0, 0d0, 0d0/
      data s/ 0d0, 1d0, 0d0, 0d0, 1d0, 0d0/
      data t/ 1d0, 1d0, 1d0,-1d0,-1d0,-1d0/
      data u/ 0d0, 0d0, 1d0, 0d0, 0d0, 1d0/
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
      integer nen,ngauss,nec,nee,i,l,nshsize,ngssize
      double precision rr,ss,tt,uu,drr,dss,dtt
      double precision trivol
c
c...  definitions
c
      trivol=half
      nshsize=(nsd+1)*nenmax*ngaussmax
      ngssize=(nsd+1)*ngaussmax
c
c...  Linear wedge definition
c
      nen=isix
      ngauss=ione
      nec=nsd*nen
      nee=ndof*nen
      gauss(1,1)=zero
      gauss(2,1)=zero
      gauss(3,1)=zero
      gauss(4,1)=two*trivol
      if(intord.ne.2) then
        ngauss=itwo
        do l=1,ngauss
          gauss(1,l)=third
          gauss(2,l)=third
          gauss(3,l)=root3i
          gauss(4,l)=trivol
        end do
        gauss(3,2)=-root3i
      end if 
c
      infetype(1)=ngauss
      infetype(2)=nen
      infetype(3)=nec
      infetype(4)=nee
c
      do l=1,ngauss
        do i=1,nen
          rr=one-r(i)+r(i)*gauss(1,l)
          ss=one-s(i)+s(i)*gauss(2,l)
          tt=one+t(i)*gauss(3,l)
          uu=one-u(i)+u(i)*gauss(1,l)
          drr=r(i)-u(i)
          dss=s(i)-u(i)
          dtt=t(i)
          sh(4,i,l)=half*rr*ss*tt*uu
          sh(1,i,l)=half*drr*ss*tt
          sh(2,i,l)=half*dss*rr*tt
          sh(3,i,l)=half*dtt*rr*ss
        end do
      end do
      call dcopy(nshsize,sh,ione,shj,ione)
c
      return
      end
c
c version
c $Id: plinwedge.f,v 1.2 2004/07/07 19:34:02 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
