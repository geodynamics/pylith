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
      subroutine read_timdat(delt,alfa,utol,ftol,etol,times,tunits,
     & maxstp,maxit,ntdinit,lgdef,itmax,nintg,lastep,kr,tfile,
     & ierr,errstrng)
c
c...program to read in time step data
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         3:  Read error
c         5:  Time units not specified
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
      integer nintg,lastep,kr,ierr
      integer maxstp(nintg),maxit(nintg),ntdinit(nintg),lgdef(nintg)
      integer itmax(nintg)
      double precision delt(nintg),alfa(nintg),utol(nintg),ftol(nintg)
      double precision etol(nintg),times(lastep+1)
      double precision tunits
      character tfile*(*),errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic index
c
c...  local variables
c
      integer i,j,n,nstep,nc
      character dummy*80
c
cdebug      write(6,*) "Hello from read_timdat_f!"
c
      ierr=izero
c
c...  open input file and skip over unit definitions
c
      open(kr,file=tfile,status="old",err=20)
      call pskip(kr)
      read(kr,"(a80)") dummy
      i=index(dummy,"=")
      if(i.eq.izero) then
        ierr=5
        errstrng="read_timdat"
        return
      end if
      call pskip(kr)
c
c...  read information on time step groups
c
      do i=1,nintg
        read(kr,*,err=30,end=30) n,maxstp(i),delt(i),alfa(i),
     &   maxit(i),ntdinit(i),lgdef(i),utol(i),ftol(i),etol(i),
     &   itmax(i)
        if(utol(i).le.zero) utol(i)=1.d-7
        if(ftol(i).le.zero) ftol(i)=1.d-4
        if(etol(i).le.zero) etol(i)=1.d-7
        delt(i)=tunits*delt(i)
      end do
      close(kr)
c
      delt(1)=zero
      maxstp(1)=izero
c
c...  create times array
c
      times(1)=zero
      nstep=izero
      do i=2,nintg
        do j=1,maxstp(i)
          nstep=nstep+ione
          nc=nstep+1
          times(nc)=times(nstep)+delt(i)
        end do
      end do
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        errstrng="read_timdat"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_timdat"
        close(kr)
        return
c
      end
c
c version
c $Id: read_timdat.f,v 1.4 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
