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
      subroutine cklock(iwink,histry,ltim,lngth,nstep,nhist,lastep,
     & unlck)
c
c...checks whether fault locking forces are being imposed or removed
c   in the current time step.  if so, rff is forced by setting
c   ltim=.true.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      double precision cutoff
      parameter(cutoff=0.8d0)
c
c...  subroutine arguments
c
      integer lngth,nstep,nhist,lastep
      integer iwink(2,lngth)
      double precision histry(nhist,lastep+1)
      logical ltim,unlck
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer k,mode,ihist
      double precision diff
      logical test
c
cdebug      write(6,*) "Hello from cklock_f!"
c
      unlck=.false.
      test=.false.
      do k=1,lngth
        mode=iwink(1,k)
        if(mode.lt.0) then
          ihist=-mode
          diff=histry(ihist,nstep+1)-histry(ihist,nstep)
          if(abs(diff).gt.cutoff) ltim=.true.
          if(abs(diff).gt.cutoff) test=.true.
          if(diff.lt.-cutoff) unlck=.true.
        end if
        if(test.and.unlck) go to 10
      end do
10    continue
      return
      end
c
c version
c $Id: cklock.f,v 1.2 2004/06/21 19:09:38 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
