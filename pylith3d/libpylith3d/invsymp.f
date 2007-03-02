c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2007 Rensselaer Polytechnic Institute
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
      subroutine invsymp(nsize,ndim,mat,ierr,errstrng)
c
c...  routine to provide the inverse of a symmetric, positive definite
c     matrix stored in packed form.
c     This routine uses LAPACK routines.
c
      implicit none
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nsize,ndim,ierr
      double precision mat(nsize)
      character errstrng*(*)
c
cdebug      write(6,*) "Hello from invsymp_f!"
c
      call dpptrf('u',ndim,mat,ierr)
      if(ierr.ne.izero) then
        errstrng="invsymp (1)"
        if(ierr.lt.izero) then
          ierr=117
        else
          ierr=118
        end if
        return
      end if
      call dpptri('u',ndim,mat,ierr)
      if(ierr.ne.izero) then
        errstrng="invsymp (2)"
        if(ierr.lt.izero) then
          ierr=117
        else
          ierr=119
        end if
      end if
      return
      end
