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
      subroutine write_tractions(tractionverts,tractionvals,
     & numtractions,nsnodes,kw,idout,ofile,ierr,errstrng)
c
c        program to print the traction BC info.
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
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numtractions,nsnodes,kw,idout,ierr
      integer tractionverts(nsnodes,numtractions)
      double precision tractionvals(ndof,numtractions)
      character ofile*(*),errstrng*(*)
c
c...  local constants
c
c
c...  intrinsic functions
c
      intrinsic index
c
c...  local variables
c
      integer i,j
      logical nonzed
c
      ierr=0
      if(numtractions.eq.izero.or.idout.eq.izero) return
c
c...  write traction BC to ascii output file, if desired
c
      call open_append(kw,ofile,"old","",ierr,errstrng,
     & "write_tractions")
      if(ierr.ne.izero) return
      write(kw,1000,err=50)
      if(nsnodes.eq.3) then
        do i=1,numtractions
          write(kw,"(3i7,3(2x,1pe15.8))")
     &    (tractionverts(j,i),j=1,nsnodes),(tractionvals(j,i),j=1,ndof)
        end do
      else if(nsnodes.eq.4) then
        do i=1,numtractions
          write(kw,"(4i7,3(2x,1pe15.8))")
     &    (tractionverts(j,i),j=1,nsnodes),(tractionvals(j,i),j=1,ndof)
        end do
      end if
      close(kw)
c
c...  normal return
c
      return
c
c...  error writing to output file
c
50    continue
        ierr=4
        errstrng="write_tractions"
        close(kw)
        return
c
 1000 format(1x,///,
     &' t r a c t i o n    b o u n d a r y    c o n d i t i o n s',///,
     &'     each set of vertices (nodes) represents the face of an',/,
     &'     element to which tractions are applied.',/,
     &'     the  tractions are specified as a stress vector in the',/,
     &'     global (x,y,z) coordinate system.',//,
     &'     vertices                                 tractions',/)
      end
c
c version
c $Id: write_tractions.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
