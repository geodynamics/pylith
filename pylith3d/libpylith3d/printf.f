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
      subroutine printf(tfault,dfault,deltp,nfault,numfn,idout,
     & idsk,kw,kp)
c
c...prints out split node displacements
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c... subroutine arguments
c
      integer numfn,idout,idsk,kw,kp
      integer nfault(3,numfn)
      double precision tfault(ndof,numfn),dfault(ndof,numfn),deltp
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer npage,n,i,j
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
cdebug      write(6,*) "Hello from printf_f!"
c
      if(numfn.eq.0) return
      npage=50
      do n=1,numfn
        if((n.eq.1.or.mod(n,npage).eq.0).and.idout.gt.1) then
          write(kw,1000) (labeld(i),i=1,ndof)
          write(kw,*) ' '
        end if
        if(idout.gt.1) write(kw,2000) nfault(1,n),nfault(2,n),
     &   (tfault(i,n),i=1,ndof)
        if(deltp.eq.zero) then
          if(idsk.eq.1) write(kp,3000) (tfault(i,n),i=1,ndof)
        else
          if(idsk.eq.1) write(kp,3000) (tfault(i,n),i=1,ndof),
     &     (dfault(j,n)/deltp,j=1,ndof)
        end if
      end do
      if(idsk.eq.2) then
        write(kp) tfault
        if(deltp.ne.zero) then
          write(kp) deltp
          write(kp) dfault
        end if
      end if
1000  format(///1x,'s p l i t   n o d e   d i s p l a c e m e n t s '//
     &    1x,'element',4x,'node',15x,6(a4,14x))
2000  format(1x,i7,4x,i7,6x,6(1pe18.5))
3000  format(12e15.7)
      return
      end
c
c version
c $Id: printf.f,v 1.3 2004/08/25 01:12:48 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
