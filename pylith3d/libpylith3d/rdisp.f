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
      subroutine rdisp(d,skew,numnp)
c
c...rotates displacement vector to global coordinate directions
c***  Note:  this could probably be done much more efficiently
c     using a skew array of dimension numrot in conjunction with
c     an index array (of the same dimension).
c     Also note:  an additional index array could keep track of which
c     skew BC correspond to slippery nodes (for the case of automatic
c     skew computation).
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numnp
      double precision d(ndof,numnp),skew(nskdim,numnp)
c
c...  local variables
c
      integer i
      double precision rot(3,3),dtemp(3)
c
      do i=1,numnp
        if((skew(1,i).ne.zero).or.(skew(nskdim,i).ne.zero)) then
          call formrt(skew(1,i),rot)
	  call dcopy(ndof,d(1,i),ione,dtemp,ione)
	  call dgemv("n",ndof,ndof,one,rot,ithree,dtemp,ione,zero,
     &     d(1,i),ione)
        end if
      end do
      return
      end
c
c version
c $Id: rdisp.f,v 1.4 2005/06/08 21:48:12 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
