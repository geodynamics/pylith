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
      subroutine formf_ss(
     & bintern,neq,                                                     ! force
     & x,numnp,iddmat,                                                  ! global
     & s,stemp,                                                         ! stiff
     & dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,                      ! elemnt
     & infmatmod,                                                       ! materl
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & skew,numrot,                                                     ! skew
     & nfault,dfault,tfault,numfn,                                      ! split
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c      generates forces due to faulted nodes
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer neq,numnp,nvfamilies,numelv,nen,ngauss,nee,numrot,numfn
      integer ierr
      integer iddmat(nstr,nstr)
      integer ien(nen,numelv),lm(ndof*nen,numelv),lmx(ndof*nen,numelv)
      integer ivfamily(5,nvfamilies),infmatmod(6,nmatmodmax)
      integer nfault(3,numfn)
      character errstrng*(*)
      double precision bintern(neq),x(nsd,numnp),s(neemax*neemax)
      double precision stemp(neemax*neemax),dfault(ndof,numfn)
      double precision tfault(ndof,numfn),dmat(nddmat*ngauss,numelv)
      double precision gauss(nsd+1,ngauss)
      double precision sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss)
      double precision skew(nskdim,numnp)
c
c...  external routines
c
      external getshape,bmatrix
c
c...  local variables
c
      integer i,ielg
      double precision p(60),dl(60)
cdebug      integer idb
c
cdebug      write(6,*) "Hello from formf_ss_f!"
c
c
c...  loop over split node entries
c
      do i=1,numfn
        ielg=nfault(1,i)
c
c...  form element stiffness matrix
c
        call formes_ss(
     &   x,numnp,iddmat,                                                ! global
     &   s,stemp,                                                       ! stiff
     &   dmat(1,ielg),ien(1,ielg),lm(1,ielg),ielg,                      ! elemnt
     &   gauss,sh,shj,nen,ngauss,nee,                                   ! eltype
     &   skew,numrot,                                                   ! skew
     &   getshape,bmatrix,                                              ! bbar
     &   ierr,errstrng)                                                 ! errcode
c
        if(ierr.ne.izero) return
c
        call lflteq(dl,dfault(1,i),nfault(1,i),ien(1,ielg),nen)
	call dsymv("u",nee,one,s,nee,dl,ione,zero,p,ione)
        call addfor(bintern,p,lm(1,ielg),lmx(1,ielg),neq,nee)
      end do
cdebug      write(6,*) "bintern:",(bintern(idb),idb=1,neq)
      return
      end
c
c version
c $Id: formf_ss.f,v 1.8 2005/04/05 23:00:22 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
