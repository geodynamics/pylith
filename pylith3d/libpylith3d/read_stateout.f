c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine read_stateout(istatout,nstatout,kr,stfile,
     & ierr,errstrng)
c
c...     reads integer array indicating which state variables
c        to output.  It is assumed that the same set of variables are
c        desired for both the elastic and time-dependent solutions.
c
c     The istatout array specifies output options for each individual
c     state variable.  At present there are a maximum of 30 possible
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
c         1:  Error opening input file
c         3:  Read error
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
      integer kr,ierr
      integer istatout(nstatesmax,3),nstatout(3)
      character stfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  local variables
c
      integer i,j
c
c...  included variable definitions
c
c
      ierr=izero
      call ifill(istatout,izero,3*nstatesmax)
      call ifill(nstatout,izero,3)
c
c...  open input file and read information on state variable output.
c     The same variables will be output for both the elastic and
c     time-dependent solutions.
c
      open(kr,file=stfile,status="old",err=20)
      do i=1,3
        call pskip(kr)
        read(kr,*,err=30,end=30) nstatout(i),
     &   (istatout(j,i),j=1,nstatout(i))
      end do
      close(kr)
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
      end
c
c version
c $Id: read_stateout.f,v 1.7 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
