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
      subroutine read_stateout(istatout,kr,kw,kp,idout,idsk,stfile,
     & ofile,pfile,ierr,errstrng)
c
c...     reads integer array indicating which state variables
c        to output.  It is assumed that the same set of variables are
c        desired for both the elastic and time-dependent solutions.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         2:  Error opening output file
c         3:  Read error
c         4:  Write error
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer kr,kw,kp,idout,idsk,ierr
      integer istatout(2,nstatesmax)
      character stfile*(*),ofile*(*),pfile*(*),errstrng*(*)
c
c...  local variables
c
      integer i,j
c
      ierr=izero
c
c...  open input file and read information on state variable output.
c     The same variables will be output for both the elastic and
c     time-dependent solutions.
c
      open(kr,file=stfile,status="old",err=20)
      call pskip(kr)
      read(kr,*,err=30,end=30) ((istatout(j,i),j=1,2),i=1,nstatesmax)
      close(kr)
c
c...  output results, if desired
c
      if(idout.gt.izero) then
        open(kw,file=ofile,err=40,status="old",access="append")
        write(kw,800,err=50) ((istatout(j,i),j=1,2),i=1,nstatesmax)
        close(kw)
      end if
      if(idsk.eq.izero) then
        open(kp,file=pfile,err=40,status="old",access="append")
        write(kp,810,err=50) ((istatout(j,i),j=1,2),i=1,nstatesmax)
      end if
      if(idsk.eq.1) then
        open(kp,file=pfile,err=40,status="old",access="append",
     &   form="unformatted")
        write(kp,err=50) istatout
      end if
      close(kp)
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        errstrng="read_stateout"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_stateout"
        close(kr)
        return
c
c...  error opening output file
c
 40   continue
        ierr=2
        errstrng="read_stateout"
        close(kw)
        close(kp)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="read_stateout"
        close(kw)
        close(kp)
        return
c
 800  format(//,
     & ' state variable output options',/,
     & ' These integers indicate which state variables are to be',/,
     & ' output for the elastic and time-dependent solutions.',/,
     & ' For stress/strain increments, either a rate (1) or an',/,
     & ' increment (2) can be selected.',//,
     & '  State Variable           Total      Increment',/,
     & '   Stress                 ',i5,8x,i5,/,
     & '   Strain                 ',i5,8x,i5,/,
     & '   Viscous strain         ',i5,8x,i5,/,
     & '   Plastic strain         ',i5,8x,i5)
 810  format(16i5)
c
      end
c
c version
c $Id: read_stateout.f,v 1.2 2004/08/24 16:45:38 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
