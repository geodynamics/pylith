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
cdebug      integer nmin,nmax
cdebug      double precision wavg
c
cdebug      write(6,*) "Hello from makemsr_f!"
c
cdebug      write(6,*) "neq,nnz,iwork,nmin,nmax,wavg:",
cdebug     & neq,nnz,iwork,nmin,nmax,wavg
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
cdebug      write(6,*) "nmin,nmax,wavg,neq:",nmin,nmax,wavg,neq
cdebug      call flush(6)
      call MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, ierr)
      call MatSeqAIJSetPreallocation(A, 0, row_array(row_offset+1),
     &                               ierr)
      call ISRestoreIndices(rowSizes, row_array, row_offset, ierr)
      call ISDestroy(rowSizes, ierr)
      return
      end
c
c version
c $Id: makemsr.F,v 1.2 2005/04/01 23:36:17 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
