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
      subroutine elas_matinit_cmp_ss(
     & A,nnz,neq,                                                       ! sparse
     & x,d,numnp,iddmat,                                                ! global
     & dx,numslp,numsn,                                                 ! slip
     & tfault,numfn,                                                    ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,               ! elemfamily
     & nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,               ! elemfamily
     & ielg,numelv,ielindx,                                             ! elemnt,orig element num
     & prop,nprop,elas_matinit,td_matinit,prestr_matinit,matchg,tminmax,! materl
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & rtimdat,ntimdat,rgiter,                                          ! timdat
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...  compute subroutine to form the d-matrix for the elastic solution
c
      include "implicit.inc"
c
c     PETSc includes
c
#include "finclude/petsc.h"
#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
      Mat A
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
      integer nnz,neq,numnp,numslp,numsn,numfn,nelfamily,nstate,nstate0
      integer nprestrflag,ipstrs,ipauto,n0states
c
c... the ielg variable below contains the global element number, while
c    the local variable ielf contains the element number within a given
c    family.
c
      integer ielg,numelv,ielindx(numelv)
      integer nprop,nen,ngauss,nee,numrot,ierr
      logical matchg
      integer iddmat(nstr,nstr)
      integer ien(nen,nelfamily),lm(ndof*nen,nelfamily)
      integer lmx(ndof*nen,nelfamily),lmf(nen,nelfamily)
      character errstrng*(*)
      double precision x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp),tfault(ndof,numfn)
      double precision s(neemax*neemax),stemp(neemax*neemax)
      double precision state(nstate,ngauss,nelfamily)
      double precision dstate(nstate,ngauss,nelfamily)
      double precision state0(nstate0,ngauss,n0states)
      double precision dmat(nddmat,ngauss,nelfamily),prop(nprop),tminmax
      double precision gauss(nsd+1,ngauss),sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss)
      double precision skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
      include "rgiter_dim.inc"
c
c...  external routines
c
      external elas_matinit,td_matinit,prestr_matinit,getshape,bmatrix
c
c...  local variables
c
      integer ielf,l
c
cdebug      integer idb
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
      include "rgiter_def.inc"
c
cdebug      write(6,*) "Hello from elas_matinit_cmp_ss_f!"
c
      tminmax=big
c
c...  loop over elements in a family and add element stiffness to global
c     stiffness.
c
      do ielf=1,nelfamily
        do l=1,ngauss
          call elas_matinit(dmat(1,l,ielf),prop,iddmat,nprop,ierr,
     &     errstrng)
cdebug          if(ielg.eq.8) then
cdebug            write(6,*) "ielf,ielg,l,dmat:",ielf,ielg,l
cdebug            write(6,*) (dmat(idb,l,ielf),idb=1,nddmat)
cdebug          end if
        end do
        call formes_ss(
     &   x,numnp,iddmat,                                                ! global
     &   s,stemp,                                                       ! stiff
     &   dmat(1,1,ielf),ien(1,ielf),lm(1,ielf),ielg,                    ! elemnt
     &   gauss,sh,shj,nen,ngauss,nee,                                   ! eltype
     &   skew,numrot,                                                   ! skew
     &   getshape,bmatrix,                                              ! bbar
     &   ierr,errstrng)                                                 ! errcod
        if(ierr.ne.izero) return
        call addstf(A,s,stemp,lm(1,ielf),lmx(1,ielf),nee,numsn,nnz,
     &   ielg,numelv,ielindx)
c       New Sieve section
cpetsc        numeq = 0
cpetsc        do i=1,nee
cpetsc          if(lm(i,ielf).ne.izero) numeq = numeq + 1
cpetsc        end do
cpetsc        k = 0
cpetsc        do i=1,nee
cpetsc          if(lm(i,ielf).ne.izero) then
cpetsc            l = 1
cpetsc            do j=1,nee
cpetsc              if(lm(j,ielf).ne.izero) then
cpetsc                stemp(k*numeq+l) = s((i-1)*nee+j)
cpetsc                l = l + 1
cpetsc              end if
cpetsc            end do
cpetsc            k = k + 1
cpetsc          end if
cpetsc        end do
cpetsc        call assembleMatrix(A, ielg, stemp, ADD_VALUES, ierr)
c
        ielg=ielg+ione
      end do
c
c     PETSc correction
c*      ielg = ielg + 1
c
      return
      end
c
c version
c $Id: elas_matinit_cmp_ss.F,v 1.5 2005/04/21 01:33:34 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
