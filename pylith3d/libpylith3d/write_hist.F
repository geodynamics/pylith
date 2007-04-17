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
      subroutine write_hist(histry,times,nhist,lastep,kw,idout,
     & ierr,errstrng)
c
c       writes load histories to the output file.
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
      integer nhist,lastep,kw,idout,ierr
      double precision histry(nhist,lastep+1),times(lastep+1)
      character errstrng*(*)
c
c...  local constants
c
      character*1 star(30)
      data star/30*'*'/
c
c...  intrinsic functions
c
      intrinsic abs,int
c
c... local variables
c
      double precision fmax,fmin,fac,range
      integer i,j,nstars,l
c
      ierr=izero
      if(idout.eq.izero) return
      write(kw,1000,err=70) nhist
      if(nhist.eq.izero) return
c
c...  write load histories
c
      do i=1,nhist
        write(kw,2000,err=70) i
c
c     echo load history and construct a miniplot
c
        fmax=-1.d32
        fmin=1.d32
        do j=1,lastep+1
          fac=histry(i,j)
          if(fac.gt.fmax) fmax=fac
          if(fac.lt.fmin) fmin=fac
        end do
        if(fmin.gt.zero) fmin=zero
        range=fmax-fmin
        do j=1,lastep+1
          if(range.ne.0) then
            nstars=int(30.0d0*((histry(i,j)-fmin)/range)+.001d0)
          else
            nstars=30
          end if
          write(kw,3000,err=70) j-1,times(j),histry(i,j),
     &     '|',(star(l),l=1,nstars)
        end do
      end do
c
c...  normal return
c
      return
c
c...  error opening output file
c
 60   continue
        ierr=2
        errstrng="write_hist"
        return
c
c...  error writing to output file
c
 70   continue
        ierr=4
        errstrng="write_hist"
        return
c
1000  format(///' l o a d   h i s t o r y   f a c t o r s'///,
     &          ' number of load histories defined (nhist) ...',i2/)
2000  format(// ' load history factor # ',i2//
     1  '  time step       time            factor',
     2 '               scaled miniplot'//)
3000  format(2x,i5,7x,1pe12.4,4x,1pe12.4,5x,31a1)
      end
c
c version
c $Id: write_hist.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
