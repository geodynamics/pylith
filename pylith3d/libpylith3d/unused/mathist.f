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
      subroutine mathist(ptmp,prop,mhist,histry,nprop,imat,nstep,nhist,
     & lastep,matchg,ierr,errstrng)
c
c...  subroutine to assign material properties based on time histories
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nprop,imat,nstep,nhist,lastep,ierr
      integer mhist(nprop)
      logical matchg
      character errstrng*(*)
      double precision ptmp(nprop),prop(nprop),histry(nhist,lastep+1)
c
c...  local variables
c
      integer i
c
cdebug      write(6,*) "Hello from mathist_f!"
c
      ierr=0
      matchg=.false.
      do i=1,nprop
        ptmp(i)=prop(i)
        if(mhist(i).ne.izero) then
          if(mhist(i).gt.nhist.or.mhist(i).lt.0) then
            ierr=100
            errstrng="mathist"
            return
          end if
          ptmp(i)=histry(mhist(i),nstep+1)*prop(i)
          matchg=.true.
        end if
      end do
c
      return
      end
c
c version
c $Id: mathist.f,v 1.1 2005/04/20 19:00:13 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
