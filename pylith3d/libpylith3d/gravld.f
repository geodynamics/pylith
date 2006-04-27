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
      subroutine gravld(p,grav,xl,iel,nen,dens,gauss,shj,ngauss,ierr,
     & errstrng)
c
c...computes the local forces due to gravitational acceleration
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
      integer iel,nen,ngauss,ierr
      character errstrng*(*)
      double precision p(ndof*nen),grav(ndof),xl(nsd,nen),dens
      double precision gauss(nsd+1,ngauss)
      double precision shj(nsd+1,nen,ngauss)
c
c...  local variables
c
      integer l,j
      double precision xs(3,3),det,rl1,rl2,rl3
c
cdebug      write(6,*) "Hello from gravld_f!"
c
      if((grav(1)*grav(1)+grav(2)*grav(2)+grav(3)*grav(3)).eq.zero)
     & return
      if(dens.eq.zero) return
      do l=1,ngauss
        call getjac(xl,xs,det,shj(1,1,l),nen,iel,ierr,errstrng)
        if(ierr.ne.izero) return
        det=det*gauss(4,l)*dens
        rl1=det*grav(1)
        rl2=det*grav(2)
        rl3=det*grav(3)
        do j=1,nen
          p(3*j-2)=p(3*j-2)+rl1*shj(4,j,l)
          p(3*j-1)=p(3*j-1)+rl2*shj(4,j,l)
          p(3*j  )=p(3*j  )+rl3*shj(4,j,l)
        end do
      end do
      return
      end
c
c version
c $Id: gravld.f,v 1.5 2005/03/21 22:13:15 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
