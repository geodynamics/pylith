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
      subroutine rstress(stn,r)
c
c...routine to rotate the stresses using the rotation matrix r
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      double precision stn(6),r(3,3)
c
c...  defined constants
c
      include "rconsts.inc"
      integer id
      parameter(id=3)
c
c...  local variables
c
      double precision stnm(3,3),tmp(3,3)
c
cdebug      write(6,*) "Hello from rstress_f!"
c
      stnm(1,1)=stn(1)
      stnm(2,2)=stn(2)
      stnm(3,3)=stn(3)
      stnm(1,2)=stn(4)
      stnm(2,3)=stn(5)
      stnm(1,3)=stn(6)
      stnm(2,1)=stnm(1,2)
      stnm(3,1)=stnm(1,3)
      stnm(3,2)=stnm(2,3)
      call dgemm("n","n",id,id,id,one,r,id,stnm,id,zero,tmp,id)
      call dgemm("n","t",id,id,id,one,tmp,id,r,id,zero,stnm,id)
      stn(1)=stnm(1,1)
      stn(2)=stnm(2,2)
      stn(3)=stnm(3,3)
      stn(4)=stnm(1,2)
      stn(5)=stnm(2,3)
      stn(6)=stnm(1,3)
      return
      end
c
c version
c $Id: rstress.f,v 1.2 2004/08/12 02:25:15 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
