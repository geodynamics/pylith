c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      SUBROUTINE zbrac(func,x1,x2,rpar,nrpar,ipar,nipar,ierr,errstrng)
c
c...  subroutine to bracket the root of a real function.  Modified
c     from Numerical Recipes.
c
c     Error codes:
c         0:  No error
c       114:  Initial bracketing values are identical
c       410:  Bracketing values not found
c
      include "implicit.inc"
c
c...  parameter definitions
c
      INTEGER NTRY
      double precision FACTOR
      PARAMETER (FACTOR=1.6d0,NTRY=50)
c
c...  subroutine arguments
c
      integer nrpar,nipar,ierr
      integer ipar(nipar)
      double precision x1,x2,rpar(nrpar)
      character errstrng*(*)
c
c...  external routines
c
      EXTERNAL func
c
c...  local variables
c
      INTEGER j
      double precision f1,df1,f2,df2
c
      ierr=0
c
      if(x1.eq.x2) then
        ierr=114
        errstrng="zbrac"
        return
      end if
      call func(x1,f1,df1,rpar,nrpar,ipar,nipar)
      call func(x2,f2,df2,rpar,nrpar,ipar,nipar)
      do j=1,NTRY
        if(f1*f2.lt.0.0d0)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+FACTOR*(x1-x2)
          call func(x1,f1,df1,rpar,nrpar,ipar,nipar)
        else
          x2=x2+FACTOR*(x2-x1)
          call func(x2,f2,df2,rpar,nrpar,ipar,nipar)
        endif
      end do
      ierr=410
      errstrng="zbrac"
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .zW.
