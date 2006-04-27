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
      integer rows, cols, localRows, localCols
      integer i,mode,k,ihist
c
cdebug      write(6,*) "Hello from winklr_f!"
c
      call MatGetSize(A, rows, cols, ierr)
      call MatGetLocalSize(A, localRows, localCols, ierr)
      call VecCreate(MPI_COMM_WORLD, winkler, ierr)
      call VecSetSizes(winkler, localRows, rows, ierr)
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
