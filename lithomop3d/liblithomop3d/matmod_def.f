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
      subroutine matmod_def(infmatmod,ismatmod)
c
c...  quick and dirty method of defining global parameters for all
c     material models.
c
c     At present, the maximum number of material models is 20, but this
c     can easily be changed by editing materials.inc.  This routine
c     defines 5 parameters for each material model:
c
c         infmatmod(1,matmod) = definition flag
c                               0 = material is not defined
c                               1 = material is defined
c         infmatmod(2,matmod) = number of state variables.  If a
c                               material is defined, it will have at
c                               least 2 state variables (stress and
c                               strain), up to a maximum of 4, if
c                               viscous or plastic strain are required.
c         infmatmod(3,matmod) = number of material properties needed to
c                               define the material.  The minimum value
c                               should be 3 for an isotropic linear
c                               elastic material (density is always the
c                               first property).
c         infmatmod(4,matmod) = material property variation flag.
c                               This flag indicates whether the material
c                               matrix varies as a function of the state
c                               variables.  If not, the same material
c                               matrix may be used for all elements in a
c                               material group.
c                               0 = material matrix independent of state
c                                   variables.
c                               1 = material matrix depends on state
c                                   variables.
c         infmatmod(5,matmod) = material property behavior flag.
c                               0 = time response of material is linear.
c                               1 = time response of material is
c                                   nonlinear.
c
c     In addition, this routine generates a simple integer (binary)
c     matrix (ismatmod) indicating whether a state variable is used for
c     a particular material model.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer infmatmod(5,nmatmodmax),ismatmod(nstatesmax,nmatmodmax)
c
c...  local variables
c
      integer i,j
c
cdebug      write(6,*) "Hello from matmod_def_f!"
c
      call ifill(infmatmod,izero,5*nmatmodmax)
      call ifill(ismatmod,izero,nstatesmax*nmatmodmax)
c
c...  Definitions for isotropic elastic material
c
      infmatmod(1,1) =ione
      infmatmod(2,1) =itwo
      infmatmod(3,1) =ithree
      ismatmod(1,1) =ione
      ismatmod(2,1) =ione
c
c...  Dummy definitions for additional elastic material models
c
      infmatmod(2,2)=itwo
      infmatmod(2,3)=itwo
      infmatmod(2,4)=itwo
      ismatmod(1,2) =ione
      ismatmod(2,2) =ione
      ismatmod(1,3) =ione
      ismatmod(2,3) =ione
      ismatmod(1,4) =ione
      ismatmod(2,4) =ione
c
c...  Dummy definition for Maxwell viscoelastic material
c
      infmatmod(2,5) = ithree
      infmatmod(3,5) = ifour
      infmatmod(4,5) = ione
      infmatmod(5,5) = izero
      ismatmod(1,5) =ione
      ismatmod(2,5) =ione
      ismatmod(3,5) =ione
c
c...  Dummy definitions for remaining materials
c
      do i=6,nmatmodmax
        do j=1,5
          infmatmod(j,i)=infmatmod(j,5)
        end do
        do j=1,nstatesmax
          ismatmod(j,i)=ismatmod(j,5)
        end do
      end do
c
      return
      end
c
c version
c $Id: matmod_def.f,v 1.2 2004/07/12 21:33:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
