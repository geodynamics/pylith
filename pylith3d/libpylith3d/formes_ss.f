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
      subroutine formes_ss(
     & x,numnp,iddmat,                                                  ! global
     & s,stemp,                                                         ! stiff
     & dmat,ien,lm,iel,                                                 ! elemnt
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...  subroutine to form the elemental stiffness matrix
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
      integer numnp,iel,nen,ngauss,nee,numrot,ierr
      integer iddmat(nstr,nstr),ien(nen),lm(ndof,nen)
      character errstrng*(*)
      double precision x(nsd,numnp),s(neemax*neemax)
      double precision stemp(neemax*neemax),dmat(nddmat,ngauss)
      double precision gauss(nsd+1,ngauss),sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss),skew(nskdim,numnp)
c
c...  external routines
c
      external getshape,bmatrix
c
c...  local variables
c
      double precision xl(60)
cdebug      integer idb
c
cdebug      write(6,*) "Hello from formes_ss_f!"
c
      call fill(s,zero,nee*nee)
c
c...  localize coordinates
c
      call lcoord(x,xl,ien,nen,numnp)
cdebug      if(iel.eq.8) then
cdebug        write(6,*) "ien:",(ien(idb),idb=1,nen)
cdebug        write(6,*) "xl:",(xl(idb),idb=1,12)
cdebug      end if
c
c...  construct local stiffness matrix, symmetrize it, and rotate for
c     skew boundary conditions
c
      call stiff_ss(
     & xl,iddmat,                                                       ! global
     & dmat,ien,iel,                                                    ! elemnt
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & s,                                                               ! stiff
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
      if(ierr.ne.izero) return
c
      if(numrot.ne.izero) call rstiff(s,stemp,skew,ien,numnp,nen,nee)
c
      return
      end
c
c version
c $Id: formes_ss.f,v 1.6 2005/04/01 23:34:13 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
