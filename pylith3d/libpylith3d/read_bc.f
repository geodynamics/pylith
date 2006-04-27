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
      subroutine read_bc(bond,dscale,vscale,fscale,ibond,numnp,
     & numbc,nconcforce,kr,bcfile,ierr,errstrng)
c
c...  subroutine to read in boundary conditions.  The bc types are
c     stored in the ibond array and the bc value is stored in the bond
c     array.
c     The input file simply consists of numbc lines with 1 + 2*ndof
c     entries per line.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         3:  Read error
c         5:  Units not specified
c       104:  BC assigned for nonexistent node
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
      integer numnp,numbc,nconcforce,kr,ierr
      integer ibond(ndof,numnp)
      double precision bond(ndof,numnp)
      double precision dscale,vscale,fscale
      character bcfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
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
      integer ihist(3),itype(3)
      integer i,j,n,imode
      double precision scale(4)
      character dummy*80
c
c...  included variable definitions
c
c
c...  open input file and define scaling factors
c
      ierr=izero
      nconcforce=izero
      open(kr,file=bcfile,status="old",err=20)
      scale(1)=zero
      scale(2)=dscale
      scale(3)=vscale
      scale(4)=fscale
c
c...  skip over units definitions.  Note that units must be specifed for
c     displacement, velocity, and force, even if that type of BC is not
c     used.
c
      do i=1,3
        call pskip(kr)
        read(kr,"(a80)") dummy
        j=index(dummy,"=")
        if(j.eq.0) then
          ierr=5
          errstrng="read_bc"
          return
        end if
      end do
c
c...  read BC, returning an error code if node is out of range
c
      call fill(bond,zero,ndof*numnp)
      call ifill(ibond,izero,ndof*numnp)
      call pskip(kr)
      do i=1,numbc
        read(kr,*,end=30,err=30) n,(ibond(j,n),j=1,ndof),
     &   (bond(j,n),j=1,ndof)
        if(n.lt.1.or.n.gt.numnp) then
          ierr=104
          errstrng="read_bc"
          return
        end if
        do j=1,ndof
          imode=ibond(j,n)
          ihist(j)=imode/10
          itype(j)=imode-10*ihist(j)
          bond(j,n)=scale(itype(j)+1)*bond(j,n)
          if(itype(j).eq.ithree) nconcforce=nconcforce+ione
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
        errstrng="read_bc"
        close(kr)
        return
c
c...  read error
c
30    continue
        ierr=3
        errstrng="read_bc"
        close(kr)
        return
c
      end
c
c version
c $Id: read_bc.f,v 1.5 2005/04/12 22:33:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
