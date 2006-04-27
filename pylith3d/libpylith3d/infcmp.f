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
      subroutine infcmp(ietypei,ietype,inf)
c
c...  routine that converts primitive element type plus infinite element
c     info to global element type.  The assumption at present is that
c     things are set up the same way as in the preshape routine.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer ietypei,ietype,inf
c
c...  local constants
c
      integer idiv(3)
      data idiv/1,10,100/
c
c...  local variables
c
      integer io(3),ind
c
cdebug      write(6,*) "Hello from infcmp_f!"
c
c
c...  simplest case:  no infinite elements
c
      if(inf.eq.izero) then
        if(ietypei.eq.ione) ietype=ione
        if(ietypei.gt.ione) ietype=ietypei+26
        if(ietypei.gt.isix) ietype=ietypei+52
      else if(ietypei.eq.ione) then
        io(3)=inf/idiv(3)
        io(2)=(inf-io(3)*idiv(3))/idiv(2)
        io(1)=inf-io(3)*idiv(3)-io(2)*idiv(2)
        ind=io(1)+ithree*io(2)+inine*io(3)
        ietype=ietypei+ind
      else if(ietypei.eq.isix) then
        io(3)=inf/idiv(3)
        io(2)=(inf-io(3)*idiv(3))/idiv(2)
        io(1)=inf-io(3)*idiv(3)-io(2)*idiv(2)
        ind=io(1)+ithree*io(2)+inine*io(3)
        ietype=ietypei+ind+26
      end if
      return
      end
c
c version
c $Id: infcmp.f,v 1.3 2005/04/08 00:41:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
