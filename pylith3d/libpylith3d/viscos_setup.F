c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
