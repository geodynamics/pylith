c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c...  Program segment to define material model 11:
c
c       Model number:                      11
c       Model name:                        ??
c       Number material properties:        ??
c       Number state variables:            ??
c       Tangent matrix varies with state:  ??
c       Material properties:               ??
c                                          ??
c                                          ??
c
      subroutine mat_prt_11(prop,nprop,matnum,idout,idsk,kw,kp,
     & ierr,errstrng)
c
c...  subroutine to output material properties for material model 11.
c
c     Error codes:
c         4:  Write error
c       101:  Attempt to use undefined material model
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nprop,matnum,idout,idsk,kw,kp,ierr
      double precision prop(nprop)
      character errstrng*(*)
c
c...  local constants
c
      character labelp(3)*15,modelname*24
      data labelp/"Density",
     &            "Young's modulus",
     &            "Poisson's ratio"/
      parameter(modelname="???????")
      integer mattype
      parameter(mattype=11)
c
c...  return error code, as this material is not yet defined
c
      ierr=101
      errstrng="mat_prt_11"
c
      return
      end
c
c
      subroutine elas_mat_11(dmat,prop,iddmat,nprop,ierr,errstrng)
c
c...  subroutine to form the material matrix for an integration point
c     for the elastic solution.  The material matrix is assumed to be
c     independent of the state variables in this case.
c     Note also that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      double precision dmat(nddmat),prop(nprop)
c
c...  return error code, as this material is not yet defined
c
      ierr=101
      errstrng="mat_prt_11"
      return
      end
c
c
      subroutine elas_strs_11(state,ee,dmat,nstate,ierr,errstrng)
c
c...  subroutine to compute stresses for the elastic solution.
c     The current total strain is contained in ee.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstate,ierr
      double precision state(nstr,nstate),ee(nstr),dmat(nddmat)
      character errstrng*(*)
c
c...  return error code, as this material is not yet defined
c
      ierr=101
      errstrng="mat_prt_11"
      return
      end
c
c
      subroutine td_matinit_11(state,dstate,dmat,prop,rtimdat,rgiter,
     & ntimdat,iddmat,tmax,nstate,nprop,matchg,ierr,errstrng)
c
c...  subroutine to form the material matrix for an integration point
c     for the time-dependent solution.  This routine is meant to be
c     called at the beginning of a time step, before strains have been
c     computed.  Thus, for some time-dependent materials, the material
c     matrix will only be an approximation of the material matrix for
c     the current iteration.
c     Note that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstate,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      double precision state(nstr,nstate),dstate(nstr,nstate)
      double precision dmat(nddmat),prop(nprop),tmax
      logical matchg
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "ntimdat_dim.inc"
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
c...  return error code, as this material is not yet defined
c
      ierr=101
      errstrng="mat_prt_11"
      return
      end
c
c
      subroutine td_strs_11(state,dstate,ee,dmat,prop,rtimdat,rgiter,
     & ntimdat,iddmat,tmax,nstate,nprop,matchg,ierr,errstrng)
c
c...  subroutine to compute the current stress for the time-dependent
c     solution.
c     Note that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstate,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      logical matchg
      double precision state(nstr,nstate),dstate(nstr,nstate),ee(nstr)
      double precision dmat(nddmat),prop(nprop),tmax
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "ntimdat_dim.inc"
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
c...  return error code, as this material is not yet defined
c
      ierr=101
      errstrng="mat_prt_11"
c
      return
      end
c
c
      subroutine td_strs_mat_11(state,dstate,ee,dmat,prop,rtimdat,
     & rgiter,ntimdat,iddmat,tmax,nstate,nprop,matchg,ierr,errstrng)
c
c...  subroutine to compute the current stress and updated material
c     matrix for the time-dependent solution.  Since this is a purely
c     elastic material, the material matrix should not change unless the
c     material properties have changed for a time step.
c     Note that only the upper triangle is used (or available), as
c     dmat is assumed to always be symmetric.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstate,nprop,ierr
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      logical matchg
      double precision state(nstr,nstate),dstate(nstr,nstate),ee(nstr)
      double precision dmat(nddmat),prop(nprop),tmax
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "ntimdat_dim.inc"
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
c...  return error code, as this material is not yet defined
c
      ierr=101
      errstrng="mat_prt_11"
c
      return
      end
c       

c version
c $Id: mat_11.f,v 1.3 2004/08/12 02:03:06 willic3 Exp $

c Generated automatically by Fortran77Mill on Tue May 18 14:18:50 2004

c End of file 
