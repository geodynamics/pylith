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
      subroutine read_diff(diforc,nslip,idhist,numslp,numdif,numnp,
     & kr,difile,ierr,errstrng)
c
c...  reads and prints differential forces applied to slippery nodes
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (if numdif.ne.zero)
c         3:  Read error
c       109:  Differential force applied to non-slippery node
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
      integer numslp,numdif,numnp,kr,ierr
      integer nslip(nsdim,numslp),idhist(numnp)
      double precision diforc(ndof,numnp)
      character difile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer n,i,j
      logical slipry
c
c...  included variable definitions
c
c
c...  open input file
c
      ierr=izero
      call fill(diforc,zero,ndof*numnp)
      call ifill(idhist,izero,numnp)
      if(numslp.eq.izero.or.numdif.eq.izero) return
      open(kr,file=difile,status="old",err=20)
c
c...  read differential force info and output results, if desired
c
      call pskip(kr)
      do i=1,numdif
        read(kr,*,end=30,err=30) n,idhist(n),(diforc(j,n),j=1,ndof)
        slipry=.false.
        do j=1,numslp
          if(n.eq.nslip(2,j)) slipry=.true.
        end do
        if(.not.slipry) then
          ierr=109
          errstrng="read_diff"
          return
        end if
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
        errstrng="read_diff"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_diff"
        close(kr)
        return
c
      end
c
c version
c $Id: read_diff.f,v 1.3 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
