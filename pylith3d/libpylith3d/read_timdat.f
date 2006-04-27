c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams
c  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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
