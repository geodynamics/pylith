c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams
c  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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
      integer rows, cols, localRows, localCols
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
      call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)
      call MatGetSize(A, rows, cols, ierr)
      call MatGetLocalSize(A, localRows, localCols, ierr)
      call VecCreate(MPI_COMM_WORLD, diag, ierr)
      call VecSetSizes(diag, localRows, rows, ierr)
      call VecSetFromOptions(diag, ierr)
      call MatGetDiagonal(A, diag, ierr)
      call VecGetArray(diag, array, offset, ierr)
      iz=0
      in=0
      ierr=0
      do i=1,localRows
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
