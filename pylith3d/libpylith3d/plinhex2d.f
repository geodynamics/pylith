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
      subroutine plinhex2d(sh,gauss,nsnodes,nsgauss,intord)
c
c... Subroutine to compute surface shape functions in natural
c    coordinates, integration points, and weights for a trilinear
c    hexahedron.
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
      integer nsnodes,nsgauss,intord
      double precision sh(nsd,nsnodes,nsgauss)
      double precision gauss(nsd,nsgauss)
c
c...  local constants
c
      double precision r(4),s(4)
      data r/-1d0, 1d0, 1d0,-1d0/
      data s/-1d0,-1d0, 1d0, 1d0/
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i,l
      double precision rr,ss
c
cdebug      integer idb
c
cdebug      write(6,*) "Hello from plinhex2d_f!"
c
c...  definitions
c
c
c...  Linear hex definition
c
      gauss(1,1)=zero
      gauss(2,1)=zero
      gauss(3,1)=four
      if(intord.ne.2) then
        do l=1,nsgauss
          gauss(1,l)=r(l)*root3i
          gauss(2,l)=s(l)*root3i
          gauss(3,l)=one
        end do
      end if
c
      do l=1,nsgauss
        do i=1,nen
          rr=one+r(i)*gauss(1,l)
          ss=one+s(i)*gauss(2,l)
          sh(3,i,l)=fourth*rr*ss
          sh(1,i,l)=fourth*r(i)*ss
          sh(2,i,l)=fourth*s(i)*rr
        end do
      end do
c
      return
      end
c
c version
c $Id: plinhex2d.f,v 1.4 2005/03/22 04:45:54 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
