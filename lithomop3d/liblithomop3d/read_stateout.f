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
c     The istatout array specifies output options for each individual
c     state variable.  At present there are a maximum of 24 possible
c     state variables, and this number may increase with the addition
c     of new material models.  For each state variable, the istatout
c     array can have either a zero (no output) or a one (output) in
c     each slot.  The options are as follows:
c
c       istatout(1,i):  Output cumulative total for current time step.
c       istatout(2,i):  Output increment over current time step.
c       istatout(3,i):  Output rate computed over current time step.
c
c      Present state variables occur in groups of 6, corresponding to
c      the number of stress/strain components, although this may change
c      in the future.  The present groups are:
c
c      1-6:    Cauchy stress
c      7-12:   Total strain
c      13-18:  Viscous strain
c      18-24:  Plastic strain
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
      integer istatout(3,nstatesmax)
      character stfile*(*),ofile*(*),pfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labels_dim.inc"
c
c...  local variables
c
      integer i,j,nstatestot,nstatesinc,nstatesrate
      integer itmp(3,nstatesmax)
c
c...  included variable definitions
c
      include "labels_def.inc"
c
      ierr=izero
      call ifill(itmp,izero,3*nstatesmax)
      nstatestot=izero
      nstatesinc=izero
      nstatesrate=izero
c
c...  open input file and read information on state variable output.
c     The same variables will be output for both the elastic and
c     time-dependent solutions.
c
      open(kr,file=stfile,status="old",err=20)
      do i=1,nstatesmax
        call pskip(kr)
        read(kr,*,err=30,end=30) (istatout(j,i),j=1,3)
      end do
      close(kr)
c
c...  dump istatout array into index array
c
      do i=1,nstatesmax
        if(istatout(1,i).ne.izero) then
          nstatestot=nstatestot+ione
          itmp(1,nstatestot)=i
        end if
        if(istatout(2,i).ne.izero) then
          nstatesinc=nstatesinc+ione
          itmp(2,nstatesinc)=i
        end if
        if(istatout(3,i).ne.izero) then
          nstatesrate=nstatesrate+ione
          itmp(3,nstatesrate)=i
        end if
      end do
c
c...  output results, if desired
c
      if(idout.gt.izero) then
        open(kw,file=ofile,err=40,status="old",access="append")
        if(nstatestot.ne.izero) then
          write(kw,700,err=50)
          write(kw,730,err=50) (labels(itmp(1,i)),i=1,nstatestot)
        end if
        if(nstatesinc.ne.izero) then
          write(kw,710,err=50)
          write(kw,730,err=50) (labels(itmp(2,i)),i=1,nstatesinc)
        end if
        if(nstatesrate.ne.izero) then
          write(kw,720,err=50)
          write(kw,730,err=50) (labels(itmp(3,i)),i=1,nstatesrate)
        end if
        close(kw)
      end if
      if(idsk.eq.ione) then
        open(kp,file=pfile,err=40,status="old",access="append")
        write(kp,810,err=50) ((istatout(j,i),j=1,3),i=1,nstatesmax)
      end if
      if(idsk.eq.itwo) then
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
 700  format(//,
     & " State variable total values to be output:",/)
 710  format(//,
     & " State variable incremental values to be output:",/)
 720  format(//,
     & " State variable rate values to be output:",/)
 730  format(4x,6(:2x,a6))
c
 810  format(16i5)
c
      end
c
c version
c $Id: read_stateout.f,v 1.4 2005/03/22 23:52:43 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
