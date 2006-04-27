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
      subroutine printd(d,deld,deltp,kout,numnp,nout,iflag,idout,
     & idsk,kto,kw,kp)
c
c...program to print displacements
c
c      options:
c
c       iflag=1, displacements are to be output
c       iflag=2, slippery node displacements are to be output using
c                entries in kout to determine which nodes to use.
c
c*  Note:  for full output, call this routine with nout = numnp.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numnp,nout,iflag,idout,idsk,kto,kw,kp
      integer kout(*)
      double precision d(ndof,numnp),deld(ndof,numnp),deltp
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
      character head(2)*53
      data head/' d i s p l a c e m e n t s                           ',
     &          ' d i f f e r e n t i a l    d i s p l a c e m e n t s'/
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer npage,n,i,m,j
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
cdebug      write(6,*) "Hello from printd_f!"
c
      if(nout.eq.0) return
      npage=50
      if(idout.gt.1.or.idsk.eq.1.or.iflag.eq.2) then
        do n=1,nout
          if((n.eq.1.or.mod(n,npage).eq.0).and.idout.gt.1) then
            write(kw,1000) head(iflag),(labeld(i),i=1,ndof)
            write(kw,*) ' '
          end if
          m=n
          if(iflag.eq.2) m=kout(n)
          if(idout.gt.1) write(kw,2000)  m,(d(i,m),i=1,ndof)
          if(deltp.eq.zero) then
            if(idsk.eq.1) write(kp,3000) m,(d(i,m),i=1,ndof)
            if(idsk.eq.2.and.iflag.eq.2) write(kp) m,(d(i,m),i=1,ndof)
          else
            if(idsk.eq.1) write(kp,3000) m,(d(i,m),i=1,ndof),
     &       (deld(j,m)/deltp,j=1,ndof)
            if(idsk.eq.2.and.iflag.eq.2) write(kp) m,
     &       (d(i,m),i=1,ndof),(deld(j,m)/deltp,j=1,ndof)
          end if
        end do
      end if
      if(idsk.eq.2.and.iflag.eq.1) then
        write(kp) d
        if(deltp.ne.zero) then
          write(kp) deltp
          write(kp) deld
        end if
      end if
      return
1000  format(///1x,a53///4x,'node number',15x,6(a4,14x))
2000  format(6x,i7,10x,6(1pe18.5))
3000  format(i7,12e15.7)
      end
c
c version
c $Id: printd.f,v 1.3 2004/08/25 01:12:48 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
