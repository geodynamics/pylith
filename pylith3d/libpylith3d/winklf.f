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
      subroutine winklf(bwink,bintern,deld,iwink,wink,histry,nwink,
     & nhist,nstep,neq,lastep)
c
c       program to compute winkler restoring forces from the
c       displacements and add them to the internal force vector bintern.
c       Note:  This method differs from the previous method since the
c       Winkler forces are now included as part of the internal force
c       vector.  Thus, the forces are now positive.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nwink,nhist,nstep,neq,lastep
      integer iwink(2,nwink)
      double precision bwink(neq),bintern(neq),deld(neq),wink(nwink)
      double precision histry(nhist,lastep+1)
c
c...  local variables
c
      integer i,mode,k,ihist
c
cdebug      write(6,*) "Hello from winklf_f!"
c
      do i=1,nwink
        mode=iwink(1,i)
        k=iwink(2,i)
        if(mode.eq.ione) then
          bwink(k)=bwink(k)+wink(i)*deld(k)
        else
          ihist=-mode
          bwink(k)=bwink(k)+wink(i)*deld(k)*histry(ihist,nstep+1)
        end if
        bintern(k)=bintern(k)+bwink(k)
      end do
      return
      end
c
c version
c $Id: winklf.f,v 1.4 2005/05/03 18:41:29 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
