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
      subroutine addstf(A,s,lm,lmx,nee,numsn,nnz)
c
c...  subroutine to add local stiffness matrix to global matrix, where
c     the global matrix is stored in modified sparse row (MSR) format.
c
      include "implicit.inc"
c
c     PETSc includes
c
#include "finclude/petsc.h"
#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
      Mat A
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nee,numsn,nnz,ierr
      integer lm(nee),lmx(nee)
      double precision s(nee,nee)
c
c...  intrinsic functions
c
      intrinsic abs,sign,dble
c
c... local variables
c
      integer nloop,i,irow,ii,j,icol,jj,itmp,l
      double precision signi,signj
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from addstf_f!"
c
cdebug      do idb=1,nee
cdebug        write(6,*) "s:",(s(jdb,idb),jdb=1,nee)
cdebug      end do
      do i=1,nee
        lm(i) = lm(i) - 1
        lmx(i) = lmx(i) - 1
      end do
      call MatSetValues(A, nee, lm, nee, lm, s, ADD_VALUES, ierr)
      do i=1,nee
        signi=sign(one,dble(lmx(i)))
        do j=1,nee
          s(i,j) = signi*s(i,j)
        end do
      end do
      call MatSetValues(A, nee, lmx, nee, lm, s, ADD_VALUES, ierr)
      do i=1,nee
        signi=sign(one,dble(lmx(i)))
        do j=1,nee
          s(i,j) = signi*s(i,j)
        end do
      end do
      do i=1,nee
        do j=1,nee
          signj=sign(one,dble(lmx(j)))
          s(i,j) = signj*s(i,j)
        end do
      end do
      call MatSetValues(A, nee, lm, nee, lmx, s, ADD_VALUES, ierr)
      do i=1,nee
        do j=1,nee
          signj=sign(one,dble(lmx(j)))
          s(i,j) = signj*s(i,j)
        end do
      end do
      do i=1,nee
        signi=sign(one,dble(lmx(i)))
        do j=1,nee
          signj=sign(one,dble(lmx(j)))
          s(i,j) = signi*signj*s(i,j)
        end do
      end do
      call MatSetValues(A, nee, lmx, nee, lmx, s, ADD_VALUES, ierr)
      do i=1,nee
        signi=sign(one,dble(lmx(i)))
        do j=1,nee
          signj=sign(one,dble(lmx(j)))
          s(i,j) = signi*signj*s(i,j)
        end do
      end do
      do i=1,nee
        lm(i) = lm(i) + 1
        lmx(i) = lmx(i) + 1
      end do
      return
      end
c
c version
c $Id: addstf.F,v 1.1 2005/03/11 04:07:40 knepley Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
