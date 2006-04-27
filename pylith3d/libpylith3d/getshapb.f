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
      subroutine getshapb(xl,sh,shj,shd,shbar,det,gauss,vol,iel,nen,
     & ngauss,ierr,errstrng)
c
c...  Subroutine to compute shape functions and derivatives at
c     Gauss points.
c
c     This is a generic routine for any element type, and it uses
c     Bbar selective integration.
c
      include "implicit.inc"
c
c...  parameter defifnitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer iel,nen,ngauss,ierr
      character errstrng*(*)
      double precision xl(nsd,nen),sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss)
      double precision shd(nsd+1,nenmax,ngaussmax),shbar(nsd+1,nenmax)
      double precision det(ngauss),gauss(nsd+1,ngauss),vol
c
c...  local variables
c
      integer l
      double precision xs(3,3)
c
c...compute shape function derivatives over number of strain integration
c   points, and then compute average dilatational component.
c
cdebug      write(6,*) "Hello from getshapb_f!"
c
      vol=zero
      do l=1,ngauss
        call getjac(xl,xs,det(l),shj(1,1,l),nen,iel,ierr,errstrng)
        if(ierr.ne.izero) return
        call getder(det(l),sh(1,1,l),shd(1,1,l),xs,nen)
        det(l)=gauss(4,l)*det(l)
        vol=vol+det(l)
      end do
      call meansh(shbar,shd,vol,det,nen,ngauss)
      return
      end
c
c version
c $Id: getshapb.f,v 1.5 2005/03/19 01:49:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
