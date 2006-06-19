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
      subroutine viscos_setup(
     & nprint,nunits,                                                   ! ioinfo
     & ofile,pfile,                                                     ! files
     & viscousStage,                                                    ! PETSc logging
     & ierr,errstrng)                                                   ! errcode
c
c...subroutine to solve the time dependent problem and perform the
c   time stepping
c
      include "implicit.inc"
c
c     PETSC includes
c
#include "finclude/petsc.h"
#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
c
c     PETSc objects
c
      PetscInt viscousStage
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer ierr
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "nprint_dim.inc"
      include "nunits_dim.inc"
c
c...  included variable definitions
c
      include "nprint_def.inc"
      include "nunits_def.inc"
c
      call PetscLogStagePush(viscousStage, ierr)
c
c...  open output files for appending, if necessary
c
      if(idout.gt.ione) then
         call open_append(kw,ofile,"old","",ierr,errstrng,"viscos")
         if(ierr.ne.izero) return
      end if
      if(idsk.eq.ione) then
         call open_append(kp,pfile,"old","",ierr,errstrng,"viscos")
         if(ierr.ne.izero) return
      else if(idsk.eq.itwo) then
         call open_append(kp,pfile,"old","unformatted",
     &        ierr,errstrng,"viscos")
         if(ierr.ne.izero) return
      end if
      return
      end
c
c version
c $Id: viscos.F,v 1.15 2005/08/05 20:17:18 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
