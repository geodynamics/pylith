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
      subroutine winklr(A,iwink,wink,histry,nstep,nwink,nhist,nnz,
     & lastep,ierr,errstrng)
c
c...program to implement winkler restoring forces on nodes
c   designated by iwink.  this program adds a constant given
c   by wink to the corresponding diagonal elements of the
c   global stiffness matrix.
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
      Vec winkler
      PetscScalar array(1)
      PetscOffset offset
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nstep,nwink,nhist,nnz,lastep,ierr
      integer iwink(2,nwink)
      double precision wink(nwink),histry(nhist,lastep+1)
      character errstrng*(*)
c
c...  local variables
c
      integer i,mode,k,ihist,rows,cols
c
cdebug      write(6,*) "Hello from winklr_f!"
c
      call MatGetSize(A, rows, cols, ierr)
      call VecCreate(MPI_COMM_WORLD, winkler, ierr)
      call VecSetSizes(winkler, PETSC_DECIDE, rows, ierr)
      call VecSetFromOptions(winkler, ierr)
      call VecGetArray(winkler, array, offset, ierr)
      do i=1,nwink
        mode=iwink(1,i)
        k=iwink(2,i)
        if(mode.eq.ione) then
          array(k+offset) = wink(i)
        else
          ihist=-mode
          if(ihist.gt.nhist) then
            ierr=100
            errstrng="winklr"
            return
          end if
          array(k+offset) = wink(i)*histry(ihist,nstep+1)
        end if
      end do
      call VecRestoreArray(winkler, array, offset, ierr)
      call MatDiagonalSet(A, winkler, ADD_VALUES, ierr)
      call VecDestroy(winkler, ierr)
      return
      end
c
c version
c $Id: winklr.F,v 1.2 2005/04/01 23:16:45 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
