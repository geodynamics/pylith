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
      subroutine get_initial_stress(
     & state,state0,infiel,nstatesz,nstatesz0,numelt,                   ! elemnt
     & infmat,infmatmod,numat,nstate0,                                  ! materl
     & infetype)                                                        ! eltype
c
c...  routine to copy computed stresses into initial stress array
c     Note:  this is an initial version.  Future versions should allow
c     arbitrary first dimensions for state and state0, and will thus
c     require material-model-specific routines to retrieve the
c     appropriate values.
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
      integer nstatesz,nstatesz0,numelt,numat,nstate0
      integer infiel(7,numelt),infmat(3,numat),infmatmod(5,nmatmodmax)
      integer infetype(4,netypes)
      double precision state(nstr,nstatesz),state0(nstate0,nstatesz0)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer iel,imat,ietype,ngauss,matmodel,nstate,indstate,l
      integer indstateg,indstate0
cdebug      integer idb,jdb
c
c...  included variable definitions
c
c
cdebug      write(6,*) "Hello from get_initial_stress_f!"
c
c
c...  loop over number of elements
c
      do iel=1,numelt
        imat=infiel(2,iel)
        ietype=infiel(3,iel)
        ngauss=infetype(1,ietype)
        matmodel=infmat(1,imat)
        nstate=infmatmod(2,matmodel)
        indstate=infiel(5,iel)
        indstate0=infiel(7,iel)
        do l=1,ngauss
          indstateg=indstate+(l-1)*nstate
          call dcopy(nstr,state(1,indstateg),ione,state0(1,indstate0),
     &     ione)
          indstate0=indstate0+ione
        end do
      end do
      return
      end
c
c version
c $Id: get_initial_stress.f,v 1.1 2005/02/23 23:47:25 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
