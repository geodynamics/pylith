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
      subroutine addsn(dl,dx,ien,lmx,nen,numnp)
c
c...adds displacements due to slip across internal interfaces to
c   the local displacement vectors
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
      integer nen,numnp
      integer ien(nen),lmx(ndof,nen)
      double precision dl(ndof,nen),dx(ndof,numnp)
c
c...  intrinsic functions
c
      intrinsic sign,dble
c
c...  local variables
c
      integer i,j,k
      double precision sgn
      logical*4 slip
c
cdebug      integer idb
c
cdebug      write(6,*) "Hello from addsn_f!"
c
      do i=1,nen
        k=ien(i)
        slip=.false.
cdebug        write(6,*) "i,k,lmx:",i,k,(lmx(idb,i),idb=1,ndof)
        do j=1,ndof
          if(lmx(j,i).ne.izero) then
            slip=.true.
            sgn=sign(one,dble(lmx(j,i)))
          end if
        end do
        if(slip) then
          do j=1,ndof
            dl(j,i)=dl(j,i)+sgn*dx(j,k)
          end do
cdebug          write(6,*) "sgn,dl:",sgn,(dl(idb,i),idb=1,ndof)
        end if
      end do
      return
      end
c
c version
c $Id: addsn.f,v 1.4 2005/04/08 00:34:37 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
