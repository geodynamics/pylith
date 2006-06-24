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
      subroutine printl(idx,iwinkx,idslp,histry,numsn,numnp,nstep,
     & nhist,nwinkx,lastep,idsk,kp)
c
c...subroutine to print array that specifies whether each slippery
c   node is locked (1) or unlocked (0) during the current time step.
c   This information is used for postprocessing.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numsn,numnp,nstep,nhist,nwinkx,lastep,idsk,kp
      integer idx(ndof,numnp),iwinkx(2,nwinkx),idslp(numsn)
      double precision histry(nhist,lastep+1)
c
c...  local variables
c
      integer i,idof,ilock,j,mode,ihist,n
c
cdebug      write(6,*) "Hello from printl_f!"
c
      if(numsn.eq.0) return
      do i=1,numsn
        n=idslp(i)
        idof=0
        ilock=1
        do j=1,ndof
          if(idx(j,n).ne.0) idof=idx(j,n)
        end do
        if(idof.eq.0) then
          ilock=1
        else if(nwinkx.ne.0) then
          do j=1,nwinkx
            if(iwinkx(2,j).eq.idof) then
              mode=iwinkx(1,j)
              if(mode.eq.0) ilock=0
              if(mode.eq.1) ilock=1
              if(mode.lt.0) then
                ihist=-mode
                ilock=0
                if(histry(ihist,nstep+1).ne.zero) ilock=1
              end if
            end if
          end do
        else
          ilock=0
        end if
        if(idsk.eq.2) write(kp) ilock
      end do
      return
      end
c
c version
c $Id: printl.f,v 1.3 2004/08/25 01:12:48 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
