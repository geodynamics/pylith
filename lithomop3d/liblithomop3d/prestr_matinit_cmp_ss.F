c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
c  Copyright 2005 Rensselaer Polytechnic Institute.
c  All worldwide rights reserved.  A license to use, copy, modify and
c  distribute this software for non-commercial research purposes only
c  is hereby granted, provided that this copyright notice and
c  accompanying disclaimer is not modified or removed from the software.
c
c  DISCLAIMER:  The software is distributed "AS IS" without any express
c  or implied warranty, including but not limited to, any implied
c  warranties of merchantability or fitness for a particular purpose
c  or any warranty of non-infringement of any current or pending patent
c  rights.  The authors of the software make no representations about
c  the suitability of this software for any particular purpose.  The
c  entire risk as to the quality and performance of the software is with
c  the user.  Should the software prove defective, the user assumes the
c  cost of all necessary servicing, repair or correction.  In
c  particular, neither Rensselaer Polytechnic Institute, nor the authors
c  of the software are liable for any indirect, special, consequential,
c  or incidental damages related to the software, to the maximum extent
c  the law permits.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine prestr_matinit_cmp_ss(
     & A,nnz,neq,                                                       ! sparse
     & x,d,numnp,iddmat,                                                ! global
     & dx,numslp,numsn,                                                 ! slip
     & tfault,numfn,                                                    ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,nelfamily,               ! elemfamily
     & nstate,nstate0,nprestrflag,ipstrs,ipauto,n0states,ielg,          ! elemfamily
     & prop,nprop,elas_matinit,td_matinit,prestr_matinit,matchg,tminmax,! materl
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & rtimdat,ntimdat,rgiter,                                          ! timdat
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...  compute subroutine to form the d-matrix for prestress
c     auto-computation
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
      integer nprestrflag,ipstrs,ipauto,n0states,ielg,nprop
      integer nen,ngauss,nee,numrot,ierr
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
      double precision gauss(nsd+1,ngauss)
      double precision sh(nsd+1,nen,ngauss)
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
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "ntimdat_def.inc"
      include "rgiter_def.inc"
c
cdebug      write(6,*) "Hello from prestr_matinit_cmp_ss_f!"
c
      tminmax=big
c
c...  loop over elements in a family and add element stiffness to global
c     stiffness.
c
      do ielf=1,nelfamily
        do l=1,ngauss
          call prestr_matinit(dmat(1,l,ielf),prop,tpois,tyoungs,iddmat,
     &     ipauto,nprop,ierr,errstrng)
        end do
c
        call formes_ss(
     &   x,numnp,iddmat,                                                ! global
     &   s,stemp,                                                       ! stiff
     &   dmat(1,1,ielf),ien(1,ielf),lm(1,ielf),ielg,                    ! elemnt
     &   gauss,sh,shj,nen,ngauss,nee,                                   ! eltype
     &   skew,numrot,                                                   ! skew
     &   getshape,bmatrix,                                              ! bbar
     &   ierr,errstrng)                                                 ! errcod
c
        if(ierr.ne.izero) return
c
        call addstf(A,s,lm(1,ielf),lmx(1,ielf),nee,numsn,nnz)
        ielg=ielg+ione
      end do
c
      return
      end
c
c version
c $Id: prestr_matinit_cmp_ss.F,v 1.3 2005/04/01 23:34:12 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
