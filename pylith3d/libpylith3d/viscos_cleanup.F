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
      subroutine viscos_cleanup(
     & ntimdat,nunits,                                                  ! timdat
     & ierr,errstrng)                                                   ! errcode
c
c     PETSC includes
c
#include "finclude/petsc.h"
#include "finclude/petsclog.h"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer ierr
      character errstrng*(*)
c
c...  included dimension and type statements
c
      include "ntimdat_dim.inc"
      include "nunits_dim.inc"
c
c...  included variable definitions
c
      include "ntimdat_def.inc"
      include "nunits_def.inc"
c
c...  close units
c
      if(idout.gt.ione) write(kw,800) ntimdat(6),ntimdat(7),ntimdat(8)
      if(idout.gt.ione) close(kw)
      close(kp)
c
 800  format(/," Total number of equilibrium iterations        = ",i7,/,
     &         " Total number of stiffness matrix reformations = ",i7,/,
     &         " Total number of displacement subiterations    = ",i7)
      call PetscLogStagePop(ierr)
      return
      end
c
c version
c $Id: viscos.F,v 1.15 2005/08/05 20:17:18 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
