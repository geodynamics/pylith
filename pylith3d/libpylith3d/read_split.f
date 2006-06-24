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
      subroutine read_split(fault,nfault,numfn,numflt,numnp,numelv,
     & kr,spfile,ierr,errstrng)
c
c...  reads and prints data on split nodes
c
c      fault(ndof,numfn)  = amount of splitting for each dof
c      nfault(3,numfn)    = element, node, and history for each entry
c      numfn              = number of split node entries
c      numflt             = total number of split nodes
c
c     Once the total number of split nodes has been computed, routine
c     id_split should be called to store the entries in idftn.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if numfn.ne.zero)
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
      integer numfn,numflt,numnp,numelv,kr,ierr
      integer nfault(3,numfn)
      double precision fault(ndof,numfn)
      character spfile*(*),errstrng*(*)
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n
c
c...  open input file and read entries
c
      ierr=izero
      numflt=izero
      if(numfn.ne.izero) then
        call fill(fault,zero,ndof*numfn)
        call ifill(nfault,izero,3*numfn)
        open(kr,file=spfile,status="old",err=20)
        call pskip(kr)
        do i=1,numfn
          read(kr,*,end=30,err=30) (nfault(j,i),j=1,3),
     &     (fault(j,i),j=1,ndof)
        end do
      end if
c
c     compute total number of split nodes
c
      if(numfn.ne.izero) then
        do n=1,numnp
          do i=1,numfn
            if(n.eq.nfault(2,i)) then
              numflt=numflt+ione
              goto 40
            end if
          end do
 40       continue
        end do
      end if
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
        errstrng="read_split"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_split"
        close(kr)
        return
c
      end
c
c version
c $Id: read_split.f,v 1.5 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
