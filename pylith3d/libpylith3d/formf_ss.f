c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
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
