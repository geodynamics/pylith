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
      subroutine scan_timdat(lastep,nintg,kr,time_units,tfile,ierr,
     & errstrng)
c
c...  subroutine to perform an initial scan of the timestep groups
c     input file to determine the number of time step groups, the total
c     number of time steps, and the units used for time.
c
c     Any error code should result in an exception since time step
c     information is required even for elastic problems.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         3:  Read error
c         5:  Units not specified
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer lastep,nintg,kr,ierr
      character time_units*(*),tfile*(*),errstrng*(*)
c
c...  local constants
c
      character def(1)*18
      data def/"time_units"/
c
c...  local variables
c
      integer nget
      integer n,maxstp,maxit,maxitc,lgdef,ibbar,itmax
      double precision delt,alfa,utol,ftol,etol
      character units(1)*80
      logical units_defined(1)
c
c...  open input file
c
cdebug      write(6,*) "Hello from scan_timdat_f!"
c
      ierr=izero
      nintg=izero
      lastep=izero
      nget=ione
      open(kr,file=tfile,status="old",err=20)
c
c...  get units, returning error 2 if they aren't found.
c
      call get_units(kr,nget,units_defined,units,def,ierr,errstrng)
      if(ierr.ne.izero) return
      time_units=units(1)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) n,maxstp,delt,alfa,maxit,maxitc,lgdef,
     &   ibbar,utol,ftol,etol,itmax
        nintg=nintg+ione
        if(nintg.gt.ione) lastep=lastep+maxstp
cdebug        write(6,*) n,maxstp,delt,alfa,maxit,maxitc,lgdef,ibbar,utol,
cdebug     &   ftol,etol,itmax,nintg,lastep
        go to 40
c
c...  normal return
c
 10   continue
        if(nintg.eq.izero) then
          ierr=3
          errstrng="scan_timdat"
        end if
        close(kr)
cdebug        write(6,*) "First error code ",ierr
        return
c
c...  error opening file
c
 20   continue
        ierr=1
        errstrng="scan_timdat"
        close(kr)
cdebug        write(6,*) "Second error code ",ierr
        return
c
c...  read error
c
 30   continue
        ierr=3
        errstrng="scan_timdat"
        close(kr)
cdebug        write(6,*) "Third error code ",ierr
        return
c
      end
c
c version
c $Id: scan_timdat.f,v 1.2 2004/07/12 20:16:53 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
