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
      subroutine traction_cmp_ss(
     & btraction,neq,                                                   ! force
     & x,d,id,numnp,                                                    ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! split
     & gauss2d,sh2d,nsnodes,nsgauss,                                    ! eltype
     & tractionverts,tractionvals,ientry,                               ! traction
     & skew,numrot,                                                     ! skew
     & ierr,errstrng)                                                   ! errcode
c
c... computation routine to add tractions to load vector.
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
      integer neq,numnp,numslp,numfn,nsnodes,nsgauss,ientry
      integer numrot,ierr
      integer id(ndof,numnp),tractionverts(nsnodes)
      double precision btraction(neq),x(nsd,numnp)
      double precision d(ndof,numnp),dx(ndof,numnp),tfault(ndof,numfn)
      double precision gauss2d(nsd,nsgauss),sh2d(nsd,nsnodes,nsgauss)
      double precision tractionvals(ndof)
      double precision skew(nskdim,numnp)
      character errstrng*(*)
c
c...  local variables
c
      integer i,j,k,l
      double precision xl(nsd,4),p(ndof,4),xs(nsd-1,nsd-1),det
c
cdebug      write(6,*) "Hello from traction_cmp_ss_f!"
c
c
c...  zero local load vector
c
      call fill(p,zero,ndof*nsnodes)
c
c...localize coordinates
c
      do i=1,nsnodes
        do j=1,nsd
          xl(j,i)=x(j,tractionverts(i))
        end do
      end do
c
c...  get contributions to local load vector
c
      do l=1,nsgauss
        call getjac2d(xl,xs,det,sh2d(1,1,l),nsnodes,ientry,ierr,
     &   errstrng)
        if(ierr.ne.izero) return
        det=det*gauss2d(3,l)
        do j=1,nsnodes
          do k=1,ndof
            p(k,j)=p(k,j)+det*sh2d(3,j,l)*tractionvals(k)
          end do
        end do
      end do
c
c...  rotate body forces if necessary
c
      if(numrot.ne.izero) call rpforc(p,skew,tractionverts,numnp,
     & nsnodes)
c
c...  add forces to global load vector
c
      do i=1,nsnodes
        do j=1,ndof
          k=id(j,tractionverts(i))
          if(k.gt.0) btraction(k)=btraction(k)+p(j,i)
        end do
      end do
      return
      end
c
c version
c $Id: traction_cmp_ss.f,v 1.2 2004/08/12 01:03:16 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
