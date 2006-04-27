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
      subroutine addstf(A,s,stemp,lm,lmx,nee,numsn,nnz,ielgpetsc)
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
      integer nee,numsn,nnz,ielgpetsc
      integer lm(nee),lmx(nee)
      double precision s(neemax*neemax),stemp(neemax*neemax)
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
      call assembleMatrix(A, ielgpetsc, stemp, ADD_VALUES, ierr)
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
