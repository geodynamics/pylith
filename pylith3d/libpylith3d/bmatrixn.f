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
      subroutine bmatrixn(b,sh,shbar,nen)
c
c...computes the linear strain-displacement matrix b(nstr,ndof*nen)
c
c      This routine does not use Hughes' B-bar formulation.
c      The shbar array is ignored.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nen
      double precision b(nstr,ndof*nenmax),sh(nsd+1,nenmax)
      double precision shbar(nsd+1,nenmax)
c
c...  local variables
c
      integer k,i
c
cdebug      write(6,*) "Hello from bmatrixn_f!"
c
      k=1
      do i=1,nen
        b(1,k  )=sh(1,i)
        b(1,k+1)=zero
        b(1,k+2)=zero
        b(2,k  )=zero
        b(2,k+1)=sh(2,i)
        b(2,k+2)=zero
        b(3,k  )=zero
        b(3,k+1)=zero
        b(3,k+2)=sh(3,i)
        b(4,k  )=sh(2,i)
        b(4,k+1)=sh(1,i)
        b(4,k+2)=zero
        b(5,k  )=zero
        b(5,k+1)=sh(3,i)
        b(5,k+2)=sh(2,i)
        b(6,k  )=sh(3,i)
        b(6,k+1)=zero
        b(6,k+2)=sh(1,i)
        k=k+3
      end do
      return
      end
c
c version
c $Id: bmatrixn.f,v 1.3 2004/08/12 01:07:17 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
