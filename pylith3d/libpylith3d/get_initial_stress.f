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
     & state,state0,ivfamily,nvfamilies,numelv,nstatesz,nstatesz0,      ! elemnt
     & infmatmod,                                                       ! materl
     & nen,ngauss)                                                      ! eltype
c
c...  routine to copy computed stresses into initial stress array
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
      integer nvfamilies,numelv,nstatesz,nstatesz0,nen,ngauss
      integer ivfamily(5,nvfamilies),infmatmod(6,nmatmodmax)
      double precision state(nstatesz),state0(nstatesz0)
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
      integer ifam,nelfamily,matmodel,indstate,indstate0,nstate,nstate0
      integer ielf,l
cdebug      integer idb,jdb
c
c...  included variable definitions
c
c
cdebug      write(6,*) "Hello from get_initial_stress_f!"
c
c
c...  loop over number of element families
c
      do ifam=1,nvfamilies
        nelfamily=ivfamily(1,ifam)
        matmodel=ivfamily(2,ifam)
        indstate=ivfamily(3,ifam)
        indstate0=ivfamily(4,ifam)
        nstate=infmatmod(2,matmodel)
        nstate0=infmatmod(6,matmodel)
c
c...  loop over elements in a family and then gauss points for each
c     element
c
        do ielf=1,nelfamily
          do l=1,ngauss
            indstate=indstate+(l-1)*nstate
            indstate0=indstate0+(l-1)*nstate0
            call dcopy(nstr,state(indstate),ione,state0(indstate0),ione)
          end do
        end do
      end do
      return
      end
c
c version
c $Id: get_initial_stress.f,v 1.2 2005/03/21 19:48:57 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
