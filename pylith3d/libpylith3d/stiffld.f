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
      subroutine stiffld(shd,stn,s,det,ngauss,nen,nee)
c
c...subroutine to compute additional stiffness contribution for
c   large deformations and add it to regular stiffness
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer ngauss,nen,nee
      double precision shd(nsd+1,nenmax,ngaussmax),stn(nstr,ngauss)
      double precision s(nee,nee),det(ngauss)
c
c...  local variables
c
      integer l,i,idim
      double precision bn(ndof*ndof,ndof*nenmax)
      double precision sb(ndof*ndof,ndof*nenmax)
      double precision stnm(ndof*ndof,ndof*ndof)
c
c...construct bn matrix, then form intermediate sb=stn*bn, and finally
c   the stiffness (bn)t*stn*bn multiplied by appropriate weight for
c   integral over element
c
cdebug      write(6,*) "Hello from stiffld_f!"
c
      idim=ndof*ndof
      call fill(stnm,zero,idim*idim)
      do l=1,ngauss
        call bnmatrix(bn,shd(1,1,l),nen)
        do i=1,3
          stnm(3*(i-1)+1,3*(i-1)+1)=stn(1,l)
          stnm(3*(i-1)+2,3*(i-1)+2)=stn(2,l)
          stnm(3*(i-1)+3,3*(i-1)+3)=stn(3,l)
          stnm(3*(i-1)+1,3*(i-1)+2)=stn(4,l)
          stnm(3*(i-1)+2,3*(i-1)+1)=stn(4,l)
          stnm(3*(i-1)+2,3*(i-1)+3)=stn(5,l)
          stnm(3*(i-1)+3,3*(i-1)+2)=stn(5,l)
          stnm(3*(i-1)+1,3*(i-1)+3)=stn(6,l)
          stnm(3*(i-1)+3,3*(i-1)+1)=stn(6,l)
        end do
        call dsymm("l","l",idim,nee,det(l),stnm,idim,bn,idim,zero,sb,
     &   idim)
        call dgemm("t","n",nee,nee,idim,one,bn,idim,sb,idim,zero,s,nee)
      end do
      return
      end
c
c version
c $Id: stiffld.f,v 1.4 2004/08/12 22:53:10 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
