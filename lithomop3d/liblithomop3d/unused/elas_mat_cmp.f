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
      subroutine elas_mat_cmp(
     & dmat,infiel,iddmat,ndmatsz,numelt,                               ! elemnt
     & prop,infmat,nprop,matgpt,elas_mat,                               ! materl
     & infetype,                                                        ! eltype
     & ierr,errstrng)                                                   ! errcode
c
c...  compute subroutine to form the d-matrix for the elastic solution
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer ndmatsz,numelt,nprop,matgpt,ierr
      integer infmat(6),infiel(6,numelt),infetype(4,netypes)
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      double precision dmat(nddmat,ndmatsz),prop(nprop)
c
c...  external routines
c
      external elas_mat
c
c...  local variables
c
      integer nmatel,imatvar,iel,inddmat0,ietype,ngauss,ng,inddmatg
      integer l,ind,inddmat
c
cdebug      write(6,*) "Hello from elas_mat_cmp_f!"
c
      nmatel=infmat(2)
      imatvar=infmat(4)
c
c...  compute d-matrix for first element in the group.  If imatvar is
c     zero, there is only one d-matrix for the group, and that is all
c     that needs to be done.  Otherwise, since this is the elastic
c     solution, the same matrix can be copied into each slot.
c
      iel=infiel(4,matgpt)
      inddmat0=infiel(6,iel)
      ietype=infiel(3,iel)
      call elas_mat(dmat(1,inddmat0),prop,iddmat,nprop,ierr,errstrng)
      if(ierr.ne.izero) return
      ngauss=infetype(1,ietype)
      ng=ngauss
      if(imatvar.eq.izero) ng=ngaussmax
      inddmatg=inddmat0
      do l=2,ng
        inddmatg=inddmatg+nddmat
        call dcopy(nddmat,dmat(1,inddmat0),ione,dmat(1,inddmatg),ione)
      end do
c
c...  loop over elements in group if there is material property
c     variation for the material type.
c
      if(imatvar.ne.0) then
        do ind=matgpt+1,matgpt+nmatel-1
          iel=infiel(4,ind)
          ietype=infiel(3,iel)
          inddmat=infiel(6,iel)
          ngauss=infetype(1,ietype)
          inddmatg=inddmat
          do l=1,ngauss
            inddmatg=inddmat+nddmat
            call dcopy(nddmat,dmat(1,inddmat0),ione,dmat(1,inddmatg),
     &       ione)
          end do
        end do
      end if
      return
      end
c
c version
c $Id: elas_mat_cmp.f,v 1.1 2004/06/25 21:41:22 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
