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
      subroutine eforce(xl,sh,shj,det,gauss,evp,p,iel,nen,ngauss,
     & getshape,bmatrix,ierr,errstrng)
c
c...this subroutine computes the effective forces at each node
c   within an element: p=(b)t*evp, where p is the local force vector,
c   b is the strain-displacement matrix, and evp is the stress
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "rconsts.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer iel,nen,ngauss,ierr
      character errstrng*(*)
      double precision xl(nsd,nen),sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss),det(ngauss)
      double precision gauss(nsd+1,ngauss),evp(nstr,ngauss)
      double precision p(ndof*nen)
      external getshape,bmatrix
c
c...  local variables
c
cdebug      integer idb,jdb,kdb
      integer l,nee
      double precision shd(4,nenmax,ngaussmax),shbar(4,20),b(6,60),vol
c
cdebug      write(6,*) "Hello from eforce_f!"
c
      nee=ndof*nen
      call getshape(xl,sh,shj,shd,shbar,det,gauss,vol,iel,nen,ngauss,
     & ierr,errstrng)
cdebug      write(6,*) "shd:"
cdebug      do idb=1,ngauss
cdebug        write(6,*) "det:",det(idb)
cdebug        do jdb=1,nen
cdebug          write(6,*) (shd(kdb,jdb,idb),kdb=1,4)
cdebug        end do
cdebug      end do
      if(ierr.ne.izero) return
c
      do l=1,ngauss
        call bmatrix(b,shd(1,1,l),shbar,nen)
cdebug        do idb=1,nee
cdebug          write(6,*) "b:",(b(jdb,idb),jdb=1,nstr)
cdebug        end do
        call dgemv("t",nstr,nee,det(l),b,nstr,evp(1,l),ione,one,p,ione)
cdebug        write(6,*) "evp:",(evp(idb,l),idb=1,nstr)
cdebug        write(6,*) "p:",(p(idb),idb=1,nee)
      end do
      return
      end
c
c version
c $Id: eforce.f,v 1.6 2005/03/19 01:49:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
