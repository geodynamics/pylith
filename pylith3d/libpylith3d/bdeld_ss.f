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
      subroutine bdeld_ss(xl,dl,sh,shj,ee,det,gauss,iel,nen,nee,ngauss,
     & getshape,bmatrix,ierr,errstrng)
c
c...  Subroutine to compute strain increments in each element.
c     Strain is evaluated at the number of gauss points specified by
c     ngauss (the locations and weights are contained in gauss).
c
c     This is a generic routine for any element type, and is meant to
c     be used only for small strain problems.
c
c     Two subroutine names are passed in as arguments.  The
c     corresponding routine is called depending on whether B-bar is
c     being used.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer iel,nen,nee,ngauss,ierr
      character errstrng*(*)
      double precision xl(nsd,nen),dl(ndof*nen)
      double precision sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss),ee(nstr,ngauss)
      double precision det(ngauss),gauss(nsd+1,ngauss)
c
c...  external routines
c
      external getshape,bmatrix
c
c...  local variables
c
      integer l
      double precision shd(4,nenmax,ngaussmax),b(6,3*nenmax)
      double precision shbar(4,nenmax),vol
c
cdebug      write(6,*) "Hello from bdeld_ss_f!"
c
c
c...compute shape function derivatives over number of strain integration
c   points, and then compute average dilatational component.
c
      call getshape(xl,sh,shj,shd,shbar,det,gauss,vol,iel,nen,ngauss,
     & ierr,errstrng)
      if(ierr.ne.izero) return
c
c...construct b matrix and compute strains at gauss point(s)
c
      do l=1,ngauss
        call bmatrix(b,shd(1,1,l),shbar,nen)
        call dgemv("n",nstr,nee,one,b,nstr,dl,ione,zero,ee(1,l),ione)
      end do
      return
      end
c
c version
c $Id: bdeld_ss.f,v 1.5 2005/03/19 01:49:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
