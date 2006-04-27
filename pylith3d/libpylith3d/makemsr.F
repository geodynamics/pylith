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
      subroutine makemsr(A,indx,link,nbrs,neq,nnz,iwork,nmin,nmax,wavg)
c
c      program to transform linked list into modified sparse row format
c
      include "implicit.inc"
c
c     PETSc includes
c
#include "finclude/petsc.h"
#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscis.h"
#include "finclude/petscmat.h"
      Mat A
      IS rowSizes
      integer row_array(1)
      PetscOffset row_offset
      PetscErrorCode ierr
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer neq,nnz,iwork,nmin,nmax
      integer indx(neq),link(iwork),nbrs(iwork)
      double precision wavg
c
c...  intrinsic functions
c
      intrinsic min,max,dble
c
c...  local variables
c
      integer i,loc,ncol
cdebug      integer idb
c
cdebug      write(6,*) "Hello from makemsr_f!"
c
      nmin=10000
      nmax=izero
      wavg=dble(neq)
      call MatSetFromOptions(A, ierr)
      call IsCreateStride(MPI_COMM_SELF, neq, 1, 1, rowSizes, ierr)
      call ISGetIndices(rowSizes, row_array, row_offset, ierr)
      do i=1,neq
        loc = indx(i)
        ncol = izero
 30     continue
        ncol = ncol+ione
        loc = link(loc)
        if(loc.gt.izero) goto 30
        row_array(i+row_offset) = ncol + 1
        nmin=min(nmin,ncol)
        nmax=max(nmax,ncol)
        wavg=wavg+dble(ncol)
      end do
      nmin=nmin+ione
      nmax=nmax+ione
      if(neq.ne.izero) wavg=wavg/dble(neq)
cFIX  call MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, ierr)
      call MatSeqAIJSetPreallocation(A, 0, row_array(row_offset+1),
     &                               ierr)
      call ISRestoreIndices(rowSizes, row_array, row_offset, ierr)
      call ISDestroy(rowSizes, ierr)
      return
      end
c
c version
c $Id: makemsr.F,v 1.3 2005/04/16 00:43:31 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
