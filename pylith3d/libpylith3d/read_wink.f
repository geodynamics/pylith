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
      subroutine read_wink(winkdef,wscal,iwinkdef,iwinkid,nwink,
     & nwinke,kr,wfile,ierr,errstrng)
c
c....program for reading and data on winkler restoring forces
c
c          winkdef(ndof,nwinke) = values of winkler restoring spring
c                             constant, force(i,j)=-wink(i,j)*disp(i,j)
c
c          iwinkdef(ndof,nwinke) = application mode:
c                               iwink = 0, no winkler forces,
c                               iwink = 1, applied throuthout computation
c                               iwink = -n, uses load history factor n
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if nwink.ne.zero)
c         3:  Read error
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
      integer nwink,nwinke,kr,ierr
      integer iwinkdef(ndof,nwinke),iwinkid(nwinke)
      double precision winkdef(ndof,nwinke),wscal(3)
      character wfile*(*),errstrng*(*)
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
c
      ierr=izero
      call ifill(iwinkid,izero,nwinke)
      call ifill(iwinkdef,izero,ndof*nwinke)
      call fill(winkdef,zero,ndof*nwinke)
c
c...  open input file
c
      if(nwink.eq.izero) return
      open(kr,file=wfile,status="old",err=20)
c
c.......read winkler force data and output results, if desired
c
      call pskip(kr)
      do i=1,nwinke
        read(kr,*,err=30,end=30) iwinkid(i),(iwinkdef(j,i),j=1,ndof),
     &   (winkdef(j,i),j=1,ndof)
        do j=1,ndof
          winkdef(j,i)=wscal(j)*winkdef(j,i)
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
        errstrng="read_wink"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_wink"
        close(kr)
        return
c
      end
c
c version
c $Id: read_wink.f,v 1.5 2005/04/16 00:45:25 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
