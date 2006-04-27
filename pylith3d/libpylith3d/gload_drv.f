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
      subroutine gload_drv(
     & bgravity,ngravflag,grav,neq,                                     ! force
     & x,d,numnp,                                                       ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! fault
     & ien,lm,lmx,lmf,ivfamily,nvfamilies,numelv,                       ! elemnt
     & prop,infmatmod,npropsz,                                          ! materl
     & gauss,shj,nen,ngauss,nee,                                        ! eltype
     & histry,rtimdat,ntimdat,nhist,lastep,gload_cmp,                   ! timdat
     & skew,numrot,                                                     ! skew
     & ierr,errstrng)                                                   ! errcode
c
c...  driver routine to add body forces due to gravity
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
      integer ngravflag,neq,numnp,numslp,numfn,nvfamilies,numelv,npropsz
      integer nen,ngauss,nee,nhist,lastep,numrot,ierr
      integer ien(nen,numelv),lm(ndof*nen,numelv),lmx(ndof*nen,numelv)
      integer lmf(nen,numelv),ivfamily(5,nvfamilies)
      integer infmatmod(6,nmatmodmax)
      character errstrng*(*)
      double precision bgravity(ngravflag*neq),grav(ndof)
      double precision x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn),prop(npropsz)
      double precision gauss(nsd+1,ngauss)
      double precision shj(nsd+1,nen,ngauss)
      double precision histry(nhist,lastep+1),skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
c
c...  external routines
c
      external gload_cmp
c
c...  local variables
c
      integer i,ielg,ifam,nelfamily,matmodel,indprop,nprop,imat
      logical matchg
      double precision gsum,dens
      double precision ptmp(100)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from gload_drv_f!"
c
      gsum=zero
      do i=1,ndof
        gsum=gsum+grav(i)*grav(i)
      end do
      if(gsum.eq.zero) return
      call fill(bgravity,zero,neq)
      ielg=ione
c
c...  loop over element families
c
      do ifam=1,nvfamilies
        nelfamily=ivfamily(1,ifam)
        matmodel=ivfamily(2,ifam)
        indprop=ivfamily(5,ifam)
        nprop=infmatmod(3,matmodel)
        imat=ifam
        matchg=.false.
        call dcopy(nprop,prop(indprop),ione,ptmp,ione)
c
        dens=ptmp(1)
c
        call gload_cmp(
     &   bgravity,grav,neq,                                             ! force
     &   x,d,numnp,                                                     ! global
     &   dx,numslp,                                                     ! slip
     &   tfault,numfn,                                                  ! fault
     &   ien(1,ielg),lm(1,ielg),lmx(1,ielg),lmf(1,ielg),nelfamily,ielg, ! elemfamily
     &   dens,matchg,                                                   ! materl
     &   gauss,shj,nen,ngauss,nee,                                      ! eltype
     &   rtimdat,ntimdat,                                               ! timdat
     &   skew,numrot,                                                   ! skew
     &   ierr,errstrng)                                                 ! errcode
c
      end do
c
      return
      end
c
c version
c $Id: gload_drv.f,v 1.6 2005/04/16 00:40:50 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
