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
      subroutine pquadtet(sh,shj,gauss,ngauss,nen,nsd,nenmax,ngaussmax,
     & intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a quadratic tetrahedron.
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
      double precision r(10),s(10),t(10),u(10)
      data r/ 1d0, 0d0, 0d0, 0d0,0.5d0,0.0d0,0.5d0,0.5d0,0.0d0,0.0d0/
      data s/ 0d0, 1d0, 0d0, 0d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0,0.0d0/
      data t/ 0d0, 0d0, 1d0, 0d0,0.0d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0/
      data u/ 0d0, 0d0, 0d0, 1d0,0.0d0,0.0d0,0.0d0,0.5d0,0.5d0,0.5d0/
c
c...  intrinsic functions
c
      intrinsic sqrt,max
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i,l,nshsize,ngssize
      double precision g1,g2,rr,ss,tt,uu,drr,dss,dtt,duu
      double precision tetvol
c
c...  definitions
c
      tetvol=sixth
      nshsize=(nsd+1)*nenmax*ngaussmax
      ngssize=(nsd+1)*ngaussmax
c
c...  Quadratic tet definition
c
      nen=10
      if(intord.eq.2) then
        ngauss=ione
        gauss(1,1)=fourth
        gauss(2,1)=fourth
        gauss(3,1)=fourth
        gauss(4,1)=tetvol
      else
        ngauss=ifour
        g1=(five-sqrt(five))/20.0d0
        g2=(five+three*sqrt(five))/20.0d0
        do l=1,ngauss
          gauss(1,l)=g1
          gauss(2,l)=g1
          gauss(3,l)=g1
          if(l.ne.ngauss) gauss(l,l)=g2
          gauss(4,l)=fourth*tetvol
        end do
      end if
      do l=1,ngauss
        do i=1,nen
          rr=r(i)*gauss(1,l)
          ss=s(i)*gauss(2,l)
          tt=t(i)*gauss(3,l)
          uu=one-rr-ss-tt
          if(i.le.4) then
            sh(4,i,l)=rr*(two*rr-one)+ss*(two*ss-one)+tt*(two*tt-one)+
     &       uu*(two*uu-one)
            sh(1,i,l)=r(i)*(four*rr-one)-u(i)*(four*uu-one)
            sh(2,i,l)=s(i)*(four*ss-one)-u(i)*(four*uu-one)
            sh(3,i,l)=t(i)*(four*tt-one)-u(i)*(four*uu-one)
          else
            sh(4,i,l)=max(rr,one)*max(ss,one)*max(tt,one)*max(uu,one)
            drr=eight*r(i)*max(ss,one)*max(tt,one)*max(uu,one)
            dss=eight*s(i)*max(rr,one)*max(tt,one)*max(uu,one)
            dtt=eight*t(i)*max(rr,one)*max(ss,one)*max(uu,one)
            duu=-eight*u(i)*max(rr,one)*max(ss,one)*max(tt,one)
            sh(1,i,l)=drr-duu
            sh(2,i,l)=dss-duu
            sh(3,i,l)=dtt-duu
          end if
        end do
      end do
      call dcopy(nshsize,sh,ione,shj,ione)
c
      return
      end
c
c version
c $Id: pquadtet.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
