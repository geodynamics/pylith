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
      subroutine write_coords(x,numnp,kw,kp,idout,idsk,
     & ierr,errstrng)
c
c...  subroutine to write nodal coordinates to plot and/or ascii files.
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
      include "ndimens.inc"
      include "nconsts.inc"
      integer npage
      parameter(npage=50)
c
c...  subroutine arguments
c
      integer numnp,kw,kp,idout,idsk,ierr
      double precision x(nsd,numnp)
      character errstrng*(*)
c
c...  included dimension and type statements
c
      include "labelc_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer i,j
c
c...  included variable definitions
c
      include "labelc_def.inc"
c
c...  output coordinates to plot file, if requested
c
      if(idsk.eq.1) then
        do i=1,numnp
          write(kp,"(i7,3(1pe16.8))",err=60) i,(x(j,i),j=1,nsd)
        end do
      else if(idsk.eq.2) then
        write(kp,err=60) x
      end if
c
c...  output coordinates to human-readable file
c
      if(idout.gt.izero) then
        do i=1,numnp
          if(i.eq.1.or.mod(i,npage).eq.0) write(kw,1000)
     &     (labelc(j),j=1,nsd)
          write(kw,"(6x,i7,10x,3(1pe20.8))",err=60) i,(x(j,i),j=1,nsd)
        end do
      end if
c
c...  normal return
c
      return
c
c...  write error
c
 60   continue
        ierr=4
        errstrng="write_coords"
        return
c
1000  format(1x,///,' n o d a l   c o o r d i n a t e   d a t a',///,
     & 4x,' node number ',13x,3(a4,18x),//)
c
      end
c
c version
c $Id: write_coords.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
