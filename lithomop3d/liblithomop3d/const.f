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
      subroutine const(maxstp,delt,alfa,maxit,maxitc,lgdef,ibbar,utol,
     & ftol,etol,itmax,nintg,i,naxstp,nfirst,rtimdat,deltp,alfap,
     & ntimdat,nstep,maxitp,maxitcp,lgdefp,ibbarp,itmaxp,gtol)
c
c...this subroutine defines the constants for each time step group.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nintg,i,naxstp,nfirst
      integer maxstp(nintg),maxit(nintg),maxitc(nintg)
      integer lgdef(nintg),ibbar(nintg),itmax(nintg)
      double precision delt(nintg),alfa(nintg),utol(nintg),ftol(nintg)
      double precision etol(nintg)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
      include "gtol_dim.inc"
c
c...  local variables
c
      integer ii
c
cdebug      write(6,*) "Hello from const_f!"
c
      ii=i
      deltp=delt(ii)
      rtimdat(1)=deltp
      alfap=alfa(ii)
      rtimdat(2)=alfap
      ntimdat(1)=nstep
      maxitp=maxit(ii)
      ntimdat(2)=maxitp
      maxitcp=maxitc(ii)
      ntimdat(3)=maxitcp
      lgdefp=lgdef(ii)
      ntimdat(4)=lgdefp
      ibbarp=ibbar(ii)
      ntimdat(5)=ibbarp
      itmaxp=itmax(ii)
      ntimdat(6)=itmaxp
      gtol(1)=utol(ii)
      gtol(2)=ftol(ii)
      gtol(3)=etol(ii)
      nfirst=naxstp+1
      naxstp=nfirst+maxstp(ii)-1
      return
      end
c
c version
c $Id: const.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
