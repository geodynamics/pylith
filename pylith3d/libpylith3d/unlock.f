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
      subroutine unlock(bwink,bintern,iwink,histry,nstep,nwink,nhist,
     & neq,lastep)
c
c       program to remove the winkler forces at a specified step.
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
      integer nstep,nwink,nhist,neq,lastep
      integer iwink(2,nwink)
      double precision bwink(neq),bintern(neq),histry(nhist,lastep+1)
c
c...  local constants
c
      double precision cutoff
      parameter(cutoff=0.8d0)
c
c...  local variables
c
      integer k,idk,mode,ihist
      double precision diff
c
cdebug      write(6,*) "Hello from unlock_f!"
c
      do k=1,nwink
        mode=iwink(1,k)
        idk=iwink(2,k)
        if(mode.lt.izero) then
          ihist=-mode
          diff=histry(ihist,nstep+1)-histry(ihist,nstep)
          if(diff.lt.-cutoff) then
            bintern(idk)=bintern(idk)-bwink(idk)
            bwink(idk)=zero
          end if
        end if
      end do
      return
      end
c
c version
c $Id: unlock.f,v 1.4 2005/05/03 18:44:14 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
