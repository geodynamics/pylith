c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
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
      subroutine formdf_ss(
     & bintern,neq,                                                     ! force
     & x,d,deld,numnp,                                                  ! global
     & s,stemp,                                                         ! stiff
     & dmat,ien,lm,lmx,infiel,iddmat,ndmatsz,numelt,nconsz,             ! elemnt
     & infmat,infmatmod,numat,                                          ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...program to compute forces due to kinematic boundary conditions
c***  Leave this routine essentially 'as-is' for now.  In the near
c***  future, there should be a more efficient method than having to
c***  loop over elements.
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
      integer neq,numnp,ndmatsz,numelt,nconsz,numat,numrot,ierr
      integer ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz)
      integer infiel(6,numelt),iddmat(nstr,nstr),infmat(3,numat)
      integer infmatmod(5,nmatmodmax),infetype(4,netypes)
      character errstrng*(*)
      double precision bintern(neq),x(nsd,numnp),d(ndof,numnp)
      double precision deld(ndof,numnp),s(neemax*neemax)
      double precision stemp(neemax*neemax),dmat(nddmat,ndmatsz)
      double precision gauss(nsd+1,ngaussmax,netypes)
      double precision sh(nsd+1,nenmax,ngaussmax,netypes)
      double precision shj(nsd+1,nenmax,ngaussmax,netypes)
      double precision skew(nskdim,numnp)
c
c...  external routines
c
      external getshape,bmatrix
c
c...  local variables
c
      integer matgpt,imat,matmodel,nmatel,imatvar,ngtest,ielg,iel,indien
      integer ietype,nen,inddmat,ngauss,nee,ngaussdim,i
      double precision p(60),dld(60)
cdebug      integer idb
c
cdebug      write(6,*) "Hello from formdf_ss_f!"
c
      matgpt=1
c
c...  loop over material groups
c
      do imat=1,numat
        matmodel=infmat(1,imat)
        nmatel=infmat(2,imat)
        imatvar=infmatmod(4,matmodel)
        ngtest=0
        if(imatvar.eq.izero) ngtest=ngaussmax
        do ielg=matgpt,matgpt+nmatel-1
          iel=infiel(4,ielg)
          indien=infiel(1,iel)
          ietype=infiel(3,iel)
          nen=infetype(2,ietype)
c
c...  localize displacement BC
c
          call ldisbc(dld,deld,ien(indien),lm(1,indien),nen,numnp)
          do i=1,ndof*nen
            if(dld(i).ne.zero) go to 100
          end do
          go to 150
 100      continue
c
c...  form element stiffness matrix
c
          inddmat=infiel(6,iel)
          ngauss=infetype(1,ietype)
          nee=infetype(4,ietype)
          ngaussdim=max(ngauss,ngtest)
c
          call formes_ss(
     &     x,numnp,                                                     ! global
     &     s,stemp,                                                     ! stiff
     &     dmat(1,inddmat),ien(indien),lm(1,indien),iddmat,iel,         ! elemnt
     &     gauss(1,1,ietype),sh(1,1,1,ietype),shj(1,1,1,ietype),        ! eltype
     &     ngauss,ngaussdim,nen,nee,                                    ! eltype
     &     skew,numrot,                                                 ! skew
     &     getshape,bmatrix,                                            ! bbar
     &     ierr,errstrng)                                               ! errcode
c
          if(ierr.ne.izero) return
c
c...  compute forces due to displacement BC and add to global vector
c
          call fill(p,zero,nee)
          call dsymv("u",nee,one,s,nee,dld,ione,zero,p,ione)
          call addfor(bintern,p,lm(1,indien),lmx(1,indien),neq,nee)
 150      continue
        end do
        matgpt=matgpt+nmatel
      end do
cdebug      write(6,*) "bintern:"
cdebug      write(6,*) (bintern(idb),idb=1,200)
      return
      end
c
c version
c $Id: formdf_ss.f,v 1.6 2005/01/05 22:10:29 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
