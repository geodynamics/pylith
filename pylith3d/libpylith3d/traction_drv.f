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
      subroutine traction_drv(
     & btraction,ntractflag,neq,                                        ! force
     & x,d,id,numnp,                                                    ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! split
     & gauss2d,sh2d,nsnodes,nsgauss,traction_cmp,                       ! eltype
     & tractionverts,tractionvals,numtractions,                         ! traction
     & skew,numrot,                                                     ! skew
     & ierr,errstrng)                                                   ! errcode
c
c... driver routine to add tractions to load vector.
c    note that time history info is presently missing.
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
      integer ntractflag,neq,numnp,numslp,numfn,nsnodes,nsgauss
      integer numtractions,numrot,ierr
      integer id(ndof,numnp),tractionverts(nsnodes,numtractions)
      double precision btraction(ntractflag*neq),x(nsd,numnp)
      double precision d(ndof,numnp),dx(ndof,numnp),tfault(ndof,numfn)
      double precision gauss2d(nsd,nsgauss),sh2d(nsd,nsnodes,nsgauss)
      double precision tractionvals(ndof,numtractions)
      double precision skew(nskdim,numnp)
      character errstrng*(*)
c
c...  external routines
c
      external traction_cmp
c
c...  local variables
c
      integer i
c
cdebug      write(6,*) "Hello from traction_drv_f!"
c
      if(ntractflag.eq.izero) return
      call fill(btraction,zero,neq)
c
c...loop over traction bc
c
      do i=1,numtractions
        call traction_cmp(
     &   btraction,neq,                                                 ! force
     &   x,d,id,numnp,                                                  ! global
     &   dx,numslp,                                                     ! slip
     &   tfault,numfn,                                                  ! split
     &   gauss2d,sh2d,nsnodes,nsgauss,                                  ! eltype
     &   tractionverts(1,i),tractionvals(1,i),i,                        ! traction
     &   skew,numrot,                                                   ! skew
     &   ierr,errstrng)                                                 ! errcode
      end do
      return
      end
c
c version
c $Id: traction_drv.f,v 1.2 2004/08/12 01:03:16 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
