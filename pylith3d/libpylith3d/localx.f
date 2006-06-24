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
      subroutine localx(idx,numnp,ien,lmx,numelv,nslip,numslp,nen)
c
c......subroutine to localize idx array and transfer sign information
c       to lmx array
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer numnp,numelv,numslp,nen
      integer idx(ndof,numnp),ien(nen,numelv),lmx(ndof,nen,numelv)
      integer nslip(nsdim,numslp)
c
c...  intrinsic functions
c
      intrinsic sign
c
c...  local variables
c
      integer i,j,k,ielg,node
c
cdebug      integer idb
c
cdebug      write(6,*) "Hello from localx_f!"
c
      call ifill(lmx,izero,ndof*nen*numelv)
      if(numslp.eq.izero) return
c
      do i=1,numslp
        ielg=nslip(1,i)
        node=nslip(2,i)
        do j=1,nen
          if(node.eq.ien(j,ielg)) then
            do k=1,ndof
              lmx(k,j,ielg)=sign(idx(k,node),nslip(2+k,i))
            end do
          end if
cdebug          write(6,*) "i,j,node,ielg,lmx:",
cdebug     &     i,j,node,ielg,(lmx(idb,j,ielg),idb=1,ndof)
        end do
      end do
      return
      end
c
c version
c $Id: localx.f,v 1.6 2005/04/08 00:31:37 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
