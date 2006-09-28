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
      subroutine read_skew(skew,runits,numrot,numnp,iskopt,kr,sfile,
     & ierr,errstrng)
c
c        program to read and print the array "skew" that specifies
c        local coordinate rotations at each node.  used for applying
c        boundary conditions in directions other than the global
c        coordinate system.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if numrot.ne.zero and
c             iskopt.eq.1)
c         3:  Read error
c         5:  Units not specified
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
      integer numrot,numnp,iskopt,kr,ierr
      double precision skew(nskdim,numnp)
      double precision runits
      character sfile*(*),errstrng*(*)
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
      integer i,j,n
      character dummy*80
c
      ierr=0
      call fill(skew,zero,nskdim*numnp)
      if(numrot.eq.izero) then
        if(iskopt.eq.itwo) numrot=numnp
        return
      end if
c
c...  open input file and skip over units definitions
c
      open(kr,file=sfile,status="old",err=20)
      call pskip(kr)
      read(kr,"(a80)") dummy
      j=index(dummy,"=")
      if(j.eq.izero) then
        ierr=5
        errstrng="read_skew"
        return
      end if
c
c...  read skew BC
c
      call pskip(kr)
      do i=1,numrot
        read(kr,*,err=30,end=30) n,(skew(j,n),j=1,nskdim)
        do j=1,nskdim
          skew(j,n)=runits*skew(j,n)
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
20    continue
        ierr=1
        errstrng="read_skew"
        close(kr)
        return
c
c...  error reading input file
c
30    continue
        ierr=3
        errstrng="read_skew"
        close(kr)
        return
c
      end
c
c version
c $Id: read_skew.f,v 1.4 2005/04/16 00:44:19 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
