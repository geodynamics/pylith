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
      subroutine gload_cmp_ss(
     & bgravity,grav,neq,                                               ! force
     & x,d,numnp,                                                       ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! fault
     & ien,lm,lmx,lmf,nelfamily,ielg,                                   ! elemfamily
     & dens,matchg,                                                     ! materl
     & gauss,shj,nen,ngauss,nee,                                        ! eltype
     & rtimdat,ntimdat,                                                 ! timdat
     & skew,numrot,                                                     ! skew
     & ierr,errstrng)                                                   ! errcode
c
c...  computation routine to add body forces due to gravity for small
c     strain problems
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
      integer neq,numnp,numslp,numfn,nelfamily,ielg,nen,ngauss,nee
      integer numrot,ierr
      integer ien(nen,nelfamily),lm(ndof*nen,nelfamily)
      integer lmx(ndof*nen,nelfamily),lmf(nen,nelfamily)
      character errstrng*(*)
      logical matchg
      double precision bgravity(neq),grav(ndof)
      double precision x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp),tfault(ndof,numfn),dens
      double precision gauss(nsd+1,ngauss)
      double precision shj(nsd+1,nen,ngauss)
      double precision skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
c
c...  local variables
c
      integer ielf
      double precision p(60),xl(60)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from gload_cmp_ss_f!"
c
c
c...  loop over elements in a family
c
      do ielf=1,nelfamily
c
c...  localize coordinates
c
        call lcoord(x,xl,ien(1,ielf),nen,numnp)
        call fill(p,zero,nee)
c
c...  compute element body force
c
        call gravld(p,grav,xl,ielg,nen,dens,gauss,shj,ngauss,ierr,
     &   errstrng)
c
c...  rotate body forces if necessary
c
        if(numrot.ne.izero) call rpforc(p,skew,ien(1,ielf),numnp,nen)
c
c...  add element forces to global vector (bgravity)
c
        call addfor(bgravity,p,lm(1,ielf),lmx(1,ielf),neq,nee)
        ielg=ielg+ione
      end do
      return
      end
c
c version
c $Id: gload_cmp_ss.f,v 1.6 2005/04/08 00:41:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
