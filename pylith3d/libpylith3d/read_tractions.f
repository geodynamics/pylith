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
      subroutine read_tractions(tractionverts,tractionvals,tscale,
     & numtractions,nsnodes,kr,tfile,ierr,errstrng)
c
c...  subroutine to read in traction BC
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         3:  Read error
c         5:  Units not specified
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer numtractions,nsnodes,kr,ierr
      integer tractionverts(nsnodes,numtractions)
      double precision tscale
      double precision tractionvals(ndof,numtractions)
      character tfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j
      character dummy*80
c
c...  included variable definitions
c
c
c...  open input file and skip past units specification.
c     Note the the current method of checking for units specification
c     is pretty sloppy, assuming that they were already specified during
c     the scan phase.
c
      if(numtractions.eq.0) return
      ierr=izero
      open(kr,file=tfile,status="old",err=20)
      call pskip(kr)
      read(kr,"(a80)") dummy
      i=index(dummy,"=")
      if(i.eq.izero) then
        ierr=ifive
        errstrng="read_tractions"
        return
      end if
      call pskip(kr)
c
c...  read the traction BC entries
c
      do i=1,numtractions
        call pskip(kr)
        read(kr,*,end=40,err=30) (tractionverts(j,i),j=1,nsnodes),
     &   (tractionvals(j,i),j=1,ndof)
        do j=1,ndof
          tractionvals(j,i)=tscale*tractionvals(j,i)
        end do
      end do
      close(kr)
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        errstrng="read_tractions"
        close(kr)
        return
c
c...  read error
c
 30   continue
        ierr=3
        errstrng="read_tractions"
        close(kr)
        return
c
c...  too few nodes
c
 40   continue
        ierr=105
        errstrng="read_tractions"
        close(kr)
        return
c
      end
c
c version
c $Id: read_tractions.f,v 1.4 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
