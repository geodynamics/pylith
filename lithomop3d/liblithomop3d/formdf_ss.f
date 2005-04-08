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
      subroutine formdf_ss(
     & bintern,neq,                                                     ! force
     & x,d,deld,numnp,iddmat,                                           ! global
     & s,stemp,                                                         ! stiff
     & dmat,ien,lm,lmx,ivfamily,nvfamilies,numelv,                      ! elemnt
     & infmatmod,                                                       ! materl
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
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
      integer neq,numnp,nvfamilies,numelv,nen,ngauss,nee,numrot,ierr
      integer iddmat(nstr,nstr)
      integer ien(nen,numelv),lm(ndof*nen,numelv),lmx(ndof*nen,numelv)
      integer ivfamily(5,nvfamilies),infmatmod(6,nmatmodmax)
      character errstrng*(*)
      double precision bintern(neq),x(nsd,numnp),d(ndof,numnp)
      double precision deld(ndof,numnp),s(neemax*neemax)
      double precision stemp(neemax*neemax),dmat(nddmat*ngauss,numelv)
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
      integer ifam,nelfamily,ielf,ielg,i
      double precision p(60),dld(60)
cdebug      integer idb
c
cdebug      write(6,*) "Hello from formdf_ss_f!"
c
      ielg=izero
c
c...  loop over element families
c
      do ifam=1,nvfamilies
        nelfamily=ivfamily(1,ifam)
        do ielf=1,nelfamily
          ielg=ielg+ione
c
c...  localize displacement BC
c
          call ldisbc(dld,deld,ien(1,ielg),lm(1,ielg),nen,numnp)
          do i=1,ndof*nen
            if(dld(i).ne.zero) go to 100
          end do
          go to 150
 100      continue
c
c...  form element stiffness matrix
c
          call formes_ss(
     &     x,numnp,iddmat,                                              ! global
     &     s,stemp,                                                     ! stiff
     &     dmat(1,ielg),ien(1,ielg),lm(1,ielg),ielg,                    ! elemnt
     &     gauss,sh,shj,nen,ngauss,nee,                                 ! eltype
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
          call addfor(bintern,p,lm(1,ielg),lmx(1,ielg),neq,nee)
 150      continue
        end do
      end do
cdebug      write(6,*) "bintern at end of formdf:",(bintern(idb),idb=1,neq)
      return
      end
c
c version
c $Id: formdf_ss.f,v 1.11 2005/04/08 00:32:19 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
