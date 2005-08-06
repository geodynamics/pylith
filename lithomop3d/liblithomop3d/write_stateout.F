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
      subroutine write_stateout(istatout,nstatout,kw,kp,idout,idsk,
     & ofile,pfile,ierr,errstrng)
c
c...     writes integer array indicating which state variables
c        to output.
c
c     The istatout array specifies output options for each individual
c     state variable.  At present there are a maximum of 24 possible
c     state variables, and this number may increase with the addition
c     of new material models.  There are three types of state variable
c     output:
c
c           1  Total accumulated values for the current time step
c           2  Incremental values from the previous step to the current
c           3  Rates computed from the previous step to the current
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
c      Three lines of input are required, corresponding to the three
c      types of state variable output.  For each line the user must
c      enter:
c      The number of state variables to output for this type (nstatout).
c        Note that the value of nstatout may be zero, in which case no
c        further output is needed for that line.
c      The state variables to output for this type (nstatout values).
c
c     Error codes:
c         0:  No error
c         2:  Error opening output file
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
      integer kw,kp,idout,idsk,ierr
      integer istatout(nstatesmax,3),nstatout(3)
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labels_dim.inc"
c
c...  local variables
c
      integer i,j
c
c...  included variable definitions
c
      include "labels_def.inc"
c
      ierr=izero
c
c...  output results, if desired
c
      if(idout.gt.izero) then
c
c...  at present, it is assumed that all f77 compilers support
c     access="append" and all f95 compilers support position="append".
c     If exceptions are found we will have to use FORTRAN_COMPILER_NAME
c     to differentiate.
c
#if defined FORTRAN_TYPE_F77
        open(kw,file=ofile,err=40,status="old",access="append")
#else
        open(kw,file=ofile,err=40,status="old",position="append")
#endif
        write(kw,700,err=50)
        if(nstatout(1).ne.izero) then
          write(kw,730,err=50) (labels(istatout(i,1)),i=1,nstatout(1))
        end if
        if(nstatout(2).ne.izero) then
          write(kw,730,err=50) 
     &     (labels(nstatesmax+istatout(i,2)),i=1,nstatout(2))
        end if
        if(nstatout(3).ne.izero) then
          write(kw,730,err=50) 
     &     (labels(2*nstatesmax+istatout(i,3)),i=1,nstatout(3))
        end if
        close(kw)
      end if
      if(idsk.eq.ione) then
#if defined FORTRAN_TYPE_F77
        open(kp,file=pfile,err=40,status="old",access="append")
#else
        open(kp,file=pfile,err=40,status="old",position="append")
#endif
        do i=1,3
          write(kp,810,err=50) nstatout(i),
     &     (istatout(j,i),j=1,nstatout(i))
        end do
        close(kp)
      end if
      if(idsk.eq.itwo) then
#if defined FORTRAN_TYPE_F77
        open(kp,file=pfile,err=40,status="old",access="append",
     &   form="unformatted")
#else
        open(kp,file=pfile,err=40,status="old",position="append",
     &   form="unformatted")
#endif
        write(kp,err=50) istatout
        close(kp)
      end if
c
c...  normal return
c
      return
c
c...  error opening output file
c
 40   continue
        ierr=2
        errstrng="write_stateout"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="write_stateout"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
        return
c
 700  format(//,
     & " State variables to be output:",/)
 730  format(4x,6(:1x,a11))
c
 810  format(16i5)
c
      end
c
c version
c $Id: write_stateout.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
