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
      subroutine write_skew(skew,numrot,iskopt,numnp,kw,idout,ofile,
     & ierr,errstrng)
c
c        program to print the array "skew" that specifies
c        local coordinate rotations at each node.  used for applying
c        boundary conditions in directions other than the global
c        coordinate system.
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
      integer numrot,iskopt,numnp,kw,idout,ierr
      double precision skew(nskdim,numnp)
      character ofile*(*),errstrng*(*)
c
c...  local constants
c
      character label(2)*5
      data label/'alpha',' beta'/
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
      if(numrot.eq.izero.or.idout.eq.izero) return
c
c...  write skew BC to ascii output file, if desired
c
      open(kw,file=ofile,err=40,status="old",access="append")
      write(kw,1000,err=50) iskopt,(label(i),i=1,nskdim)
      write(kw,1050,err=50)
      do i=1,numnp
        nonzed=.false.
        do j=1,nskdim
          if(skew(j,i).ne.zero) nonzed=.true.
        end do
        if(nonzed) write(kw,2000,err=40) i,(skew(j,i),j=1,nskdim)
      end do
      close(kw)
c
c...  normal return
c
      return
c
c...  error opening output file
c
40    continue
        ierr=2
        errstrng="write_skew"
        close(kw)
        return
c
c...  error writing to output file
c
50    continue
        ierr=4
        errstrng="write_skew"
        close(kw)
        return
c
 1000 format(1x,///,
     &' l o c a l    c o o r d i n a t e    r o t a t i o n s',///,
     &' Coordinate rotation option (iskopt) . . . . . .',i5,/,
     &'     1 = rotations are assigned',/,
     &'     2 = rotations are computed internally',///,
     &'     alpha = cc rotation (radians) from x axis in xy plane',/,
     &'      beta = cc rotation (radians) from x axis in xz plane',//,
     &'   alpha rotation is applied first, then the beta rotation--',/,
     &'   this rotates the skew coordinates into the global system.',//,
     & 1x,'node #',2(10x,a5))
 1050 format(/)
 2000 format(2x,i7,2(f15.2))
      end
c
c version
c $Id: write_skew.f,v 1.1 2005/08/05 19:58:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
