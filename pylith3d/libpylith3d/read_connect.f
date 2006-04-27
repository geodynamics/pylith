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
      subroutine read_connect(ien,mat,nen,numelv,numnp,nvfamilies,
     & kr,ifile,ierr,errstrng)
c
c      this subroutine reads element types and connectivities, material
c      types, and infinite element info.
c
c      Error codes:
c          0:  No error
c          1:  Error opening input file
c          3:  Read error
c        106:  Illegal element type
c        107:  Illegal material type
c        108:  Undefined node used for connectivity
c
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nen,numelv,numnp,nvfamilies,kr,ierr
      integer ien(nen,numelv),mat(numelv)
      character ifile*(*),errstrng*(*)
c
c...  local variables
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n,ietypev,infin
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from read_connect_f!"
c
      call ifill(ien,izero,nen*numelv)
      ierr=izero
c
c...  read connectivity, material number, and infinite element info
c
      open(kr,file=ifile,status="old",err=20)
      call pskip(kr)
      do i=1,numelv
        read(kr,*,end=30,err=30) n,ietypev,mat(i),infin,
     &   (ien(j,i),j=1,nen)
c
c...  check for illegal element type
c
        if(ietypev.le.izero.or.ietypev.gt.netypesi) then
          ierr=106
          errstrng="read_connect"
          return
        end if
clater        call infcmp(ietypev,infiel(3,i),inf)
c
c...  check for illegal material type or material model
c
        if(mat(i).le.izero.or.mat(i).gt.nvfamilies) then
          ierr=107
          errstrng="read_connect"
          return
        end if
      end do
      close(kr)
c
c......test for zero or out-of-bound entries in ien array
c
      do n=1,numelv
        do j=1,nen
          if(ien(j,n).le.izero.or.ien(j,n).gt.numnp) then
            ierr=108
            errstrng="read_connect"
            return
          end if
        end do
      end do
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        errstrng="read_connect"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_connect"
        close(kr)
        return
c
      end
c
c version
c $Id: read_connect.f,v 1.11 2005/04/13 00:37:51 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
