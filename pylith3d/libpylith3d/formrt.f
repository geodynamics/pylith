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
      subroutine formrt(skew,rot)
c
c...constructs a three dimensional rotation matrix rot from the euler
c   angles contained in the array skew.
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
      double precision skew(nskdim),rot(3,3)
c
c...  intrinsic functions
c
      intrinsic sin,cos
c
c...  local variables
c
      double precision sina,cosa,sinb,cosb
c
cdebug      write(6,*) "Hello from formrt_f!"
c
      sina=sin(skew(1))
      cosa=cos(skew(1))
      sinb=sin(skew(2))
      cosb=cos(skew(2))
      rot(1,1)=cosa*cosb
      rot(2,1)=sina*cosb
      rot(3,1)=-sinb
      rot(1,2)=-sina
      rot(2,2)=cosa
      rot(3,2)=zero
      rot(1,3)=cosa*sinb
      rot(2,3)=sina*sinb
      rot(3,3)=cosb
      return
      end
c
c version
c $Id: formrt.f,v 1.3 2004/08/12 01:22:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
