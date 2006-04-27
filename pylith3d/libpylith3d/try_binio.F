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
      subroutine try_binio(u,ierr,errstrng)
c
c...  subroutine to try binary I/O to see if it works
c
#include "binio.h"
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer u,ierr
      character errstrng*(*)
c
c...  local variables
c
      integer intlen,floatlen,p
      integer testint
      real testfloat
c
      ierr=0
c
      intlen=0
      floatlen=0
      p=0
      testint=42
      testfloat=4.2
c
#if defined FORTRAN_NO_BINARY
      ierr=400
      goto 10
#else
c
c...  open a scratch file
c
      open(u,status="scratch",
     & BINARY_OPEN_SPEC,
     & form="unformatted",err=402)
c
c...  write an integer and a float and measure their size
c
      write(u,BINARY_POS=1,BINARY_ERR(404)) testint
      inquire(u,BINARY_INQ_POS=p)
      intlen = p - 1
      write(u,BINARY_POS=p,BINARY_ERR(404)) testfloat
      inquire(u,BINARY_INQ_POS=p)
      floatlen = p - intlen - 1
      if(intlen.ne.ifour.or.floatlen.ne.ifour) then
        ierr=405
        goto 10
      end if
c
c...  read them back and verify that the values are the same
c
      testint=0
      testfloat=0.0
      read(u,BINARY_POS=1,BINARY_ERR(403)) testint
      p = intlen + 1
      read(u,BINARY_POS=p,BINARY_ERR(403)) testfloat
      close(u)
      if(testint.ne.42.or.testfloat.ne.4.2) then
        ierr=406
        goto 10
      end if
#endif
c
c...  normal return
c
      return
c
c...  error
c
 402  ierr=402
      goto 10
 403  ierr=403
      goto 10
 404  ierr=404
      goto 10
c
 10   errstrng="try_binio"
      close(u)
      return
c
      end
c
c End of file 
