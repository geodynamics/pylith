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
      subroutine write_timdat(delt,alfa,utol,ftol,etol,times,maxstp,
     & maxit,ntdinit,lgdef,itmax,nintg,lastep,kw,idout,
     & ierr,errstrng)
c
c...program to write time step data
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
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nintg,lastep,kw,idout,ierr
      integer maxstp(nintg),maxit(nintg),ntdinit(nintg),lgdef(nintg)
      integer itmax(nintg)
      double precision delt(nintg),alfa(nintg),utol(nintg),ftol(nintg)
      double precision etol(nintg),times(lastep+1)
      character errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic index
c
c...  local variables
c
      integer i,j,ii,nstep,nc
c
cdebug      write(6,*) "Hello from write_timdat_f!"
c
      ierr=izero
      if(idout.eq.izero) return
c
c...  echo input to file
c
      write(kw,1000,err=50)
      do i=1,nintg
          write(kw,2000,err=50) i,maxstp(i),delt(i),alfa(i),maxit(i),
     &   ntdinit(i),lgdef(i),utol(i),ftol(i),etol(i),itmax(i)
      end do
c
c...write out time-step/time relationship for one cycle
c
      nstep=izero
      write(kw,3000,err=50) nstep,1,times(1)
      do i=2,nintg
        ii=i
        write(kw,5000,err=50)
        do j=1,maxstp(i)
          nstep=nstep+ione
          nc=nstep+ione
          write(kw,4000,err=50) nstep,ii,times(nc)
        end do
      end do
c
c...  normal return
c
      return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="write_timdat"
        return
c
1000  format(///,' t i m e   s t e p   i n f o r m a t i o n',//,
     & '   Note:  Time step group #1 is the elastic solution',//,
     & 'group   #      step    alfa  maxit ntdinit lgdef utol',
     & '  ftol  etol  itmax',/,
     & '  #   steps    size',/)
2000  format(i3,1x,i5,2x,1pe10.4,2x,0pf4.2,2x,i5,1x,i5,1x,i5,1x,
     & 1pe7.1,1x,1pe7.1,1x,1pe7.1,1x,i5)
3000  format(//,' time-step/time correspondence:',//,
     1 ' time step #   in   group #            time',//,
     2       3x,i5,11x,i5,7x,1pe15.4,5x,'(Elastic)')
4000  format(3x,i5,11x,i5,7x,1pe15.4)
5000  format(' ')
      end
c
c version
c $Id: write_timdat.F,v 1.1 2005/08/05 19:58:07 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
