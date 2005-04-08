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
      subroutine ckdiag(A,ierr,errstrng)
c
c...program to test for zero or negative diagonals of the
c   stiffness matrix
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
      Mat A
      Vec diag
      integer rows, cols
      PetscScalar array(1)
      PetscOffset offset
c
c...  parameter definitions
c
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer ierr
      character errstrng*(*)
c
c...  local variables
c
      integer iz,in,i
c
cdebug      write(6,*) "Hello from ckdiag_f!"
c
cdebug      call MatView(A, PETSC_VIEWER_STDOUT_SELF, ierr)
      call MatGetSize(A, rows, cols, ierr)
      call VecCreate(MPI_COMM_WORLD, diag, ierr)
      call VecSetSizes(diag, PETSC_DECIDE, rows, ierr)
      call VecSetFromOptions(diag, ierr)
      call MatGetDiagonal(A, diag, ierr)
      call VecGetArray(diag, array, offset, ierr)
      iz=0
      in=0
      ierr=0
      do i=1,rows
        if(array(i+offset).eq.zero) iz=iz+1
        if(array(i+offset).lt.zero) in=in+1
cdebug        write(6,*) "i,array(i):",i,array(i+offset)
      end do
      call VecRestoreArray(diag, array, offset, ierr)
      if((iz.ne.0).or.(in.ne.0)) then
        ierr=103
        write(errstrng,1000) iz,in
      end if
 1000 format("ckdiag: ",i7," zero diags,",i7," neg diags")
      return
      end
c
c version
c $Id: ckdiag.F,v 1.3 2005/04/08 00:33:11 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
