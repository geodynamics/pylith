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
      subroutine addstf(A,s,stemp,lm,lmx,nee,numsn,nnz,ielg,
     &     numelv,ielindx)
c
c...  subroutine to add local stiffness matrix to global matrix, where
c     the global matrix is stored in modified sparse row (MSR) format.
c     This is a new version of the code that uses the assembleMatrix
c     routine from Sieve/PETSc.
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
      include "nshape.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nee,numsn,nnz
      integer lm(nee),lmx(nee)
      double precision s(neemax*neemax),stemp(neemax*neemax)
      integer ielg,numelv,ielindx(numelv),ielgnew
c
c...  intrinsic functions
c
      intrinsic abs,sign,dble
c
c... local variables
c
      integer i,j,k,l,numeq,ierr
cpetsc      integer lmp(neemax),lmxpa(neemax)
      double precision signi,signj
      logical slip
c
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from addstf_f!"
c
      slip=.false.
      do i=1,nee
cdebug        write(6,*) "i,lm(i),s(i,i):",i,lm(i),s(i,i)
        if(lm(i).ne.izero) numeq = numeq + 1
        if(lmx(i).ne.izero) slip=.true.
cpetsc        lmp(i) = lm(i) - 1
cpetsc        lmxpa(i) = abs(lmx(i))-1
      end do
c     New Sieve section
      numeq = 0
      do i=1,nee
        if(lm(i).ne.izero) numeq = numeq + 1
      end do
      k = 0
      do i=1,nee
        if(lm(i).ne.izero) then
          l = 1
          do j=1,nee
            if(lm(j).ne.izero) then
              stemp(k*numeq+l) = s((i-1)*nee+j)
              l = l + 1
            end if
          end do
          k = k + 1
        end if
      end do
      ielgnew = ielindx(ielg) - 1
      call assembleMatrix(A, ielgnew, stemp, ADD_VALUES, ierr)
cpetsc      call MatSetValues(A, nee, lmp, nee, lmp, s, ADD_VALUES, ierr)
cpetsc      if(slip) then
cpetsc        do i=1,nee
cpetsc          signi=sign(one,dble(lmx(i)))
cpetsc          do j=1,nee
cpetsc            stemp(i,j) = signi*s(i,j)
cpetsc          end do
cpetsc        end do
cpetsc        call MatSetValues(A, nee,lmp, nee,lmxpa,stemp, ADD_VALUES, ierr)
cpetsc        do i=1,nee
cpetsc          do j=1,nee
cpetsc            signj=sign(one,dble(lmx(j)))
cpetsc            stemp(i,j) = signj*s(i,j)
cpetsc          end do
cpetsc        end do
cpetsc        call MatSetValues(A, nee,lmxpa, nee,lmp,stemp, ADD_VALUES, ierr)
cpetsc        do i=1,nee
cpetsc          signi=sign(one,dble(lmx(i)))
cpetsc          do j=1,nee
cpetsc            signj=sign(one,dble(lmx(j)))
cpetsc            stemp(i,j) = signi*signj*s(i,j)
cpetsc          end do
cpetsc        end do
cpetsc        call MatSetValues(A,nee,lmxpa,nee,lmxpa,stemp, ADD_VALUES, ierr)
cpetsc      end if
      return
      end
c
c version
c $Id: addstf.F,v 1.6 2005/04/08 00:40:28 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
