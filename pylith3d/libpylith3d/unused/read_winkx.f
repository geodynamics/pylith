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
      subroutine read_winkx(winkx,wxscal,iwinkx,idx,numnp,nwinkx,
     & nwinkxe,kr,wxfile,ierr,errstrng)
c
c...  program for reading and printing data on differential winkler
c     restoring forces
c
c          winkx(nwinkx) = values of winkler restoring spring
c                             constant, force(idx(i))=-winkx(i)*disp(i)
c
c          iwtmpx(2,nwinkx) = application mode:
c                             iwinkx = 0, no winkler forces,
c                             iwinkx = 1, applied throughout computation
c                             iwinkx = -n, uses load history factor n
c
c        After input, the iwtmp array is transformed into the iwink
c        array, which has the following components:
c          iwinkx(1,i)       = nonzero application mode from iwtmp
c          iwinkx(2,i)       = global equation number to which force is
c                             applied

c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if nwinkx.ne.zero)
c         3:  Read error
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
      integer numnp,nwinkx,nwinkxe,kr,ierr
      integer iwinkx(2,nwinkx),idx(ndof,numnp)
      double precision winkx(nwinkx),wxscal(3)
      character wxfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n,nwxtot,nnz
      integer iwxtmp(3)
      double precision wxtmp(3)
c
c...  included variable definitions
c
c
c...  open output file
c
      ierr=izero
c
c...  open input file
c
      if(nwinkx.eq.izero) return
      open(kr,file=wxfile,status="old",err=20)
c
c.......read winkler force data and output results, if desired
c
      call pskip(kr)
      nwxtot=izero
      do i=1,nwinkxe
        read(kr,*,err=30,end=30) n,(iwxtmp(j),j=1,ndof),
     &   (wxtmp(j),j=1,ndof)
        nnz=izero
        do j=1,ndof
          wxtmp(j)=wxscal(j)*wxtmp(j)
          if(iwxtmp(j).ne.izero) then
            nnz=nnz+1
            nwxtot=nwxtot+1
            iwinkx(1,nwxtot)=iwxtmp(j)
            iwinkx(2,nwxtot)=idx(j,n)
            winkx(nwxtot)=wxtmp(j)
          end if
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
        errstrng="read_winkx"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_winkx"
        close(kr)
        return
c
      end
c
c version
c $Id: read_winkx.f,v 1.1 2005/04/20 19:00:14 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
