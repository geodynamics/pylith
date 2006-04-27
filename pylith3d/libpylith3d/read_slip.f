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
      subroutine read_slip(nslip,numslp,numsn,numnp,iskopt,kr,slfile,
     & ierr,errstrng)
c
c...  reads data on free slip interfaces
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if numslp.ne.zero)
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
      integer numslp,numsn,numnp,iskopt,kr,ierr
      integer nslip(nsdim,numslp)
      character slfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n
c
c...  included variable definitions
c
c
c...  open input file
c
      ierr=0
      numsn=0
      if(numslp.eq.izero) return
      call ifill(nslip,izero,nsdim*numslp)
      open(kr,file=slfile,status="old",err=20)
      call pskip(kr)
      do i=1,numslp
        read(kr,*,end=30,err=30) (nslip(j,i),j=1,nsdim)
        if(iskopt.ne.ione) then
          if(nslip(3,i).ne.izero) then
            nslip(4,i)=nslip(3,i)
            nslip(5,i)=izero
          else if(nslip(4,i).ne.izero) then
            nslip(3,i)=nslip(4,i)
            nslip(5,i)=izero
          else if(nslip(5,i).ne.izero) then
            nslip(3,i)=nslip(5,i)
            nslip(4,i)=nslip(5,i)
            nslip(5,i)=izero
          end if
        end if
      end do
c
c...  determine number of slippery nodes
c
      do n=1,numnp
        do j=1,numslp
          if(n.eq.nslip(2,j)) then
            numsn=numsn+ione
            go to 40
          end if
        end do
 40     continue
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
        errstrng="read_slip"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_slip"
        close(kr)
        return
c
      end
c
c version
c $Id: read_slip.f,v 1.6 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
