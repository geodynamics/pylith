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
      subroutine rstiff(s,stemp,skew,ien,numnp,nen,nee)
c
c...rotates the local stiffness matrix's coordinate system
c   with respect to the global coordinates at specified nodes
c   for imposition of skew boundary conditions.  k'= (r-1) k r
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numnp,nen,nee
      integer ien(nen)
      double precision s(nee,nee),stemp(nee,nee),skew(nskdim,numnp)
c
c...  local variables
c
      integer i,j,k,l,ll
      double precision rot(3,3),sum
c
c*******  New method to be implemented:
c  Use a local array (dimension 60*60?) to store BLAS results.  First loop over
c  k r multiplications, storing results in srot(1,1+ndof*(i-1)).  At the same
c  time, set values for nrot (# rotations) and irot (index of each rotation).
c  Then loop over rotated nodes.  For each rotated node, first set appropriate
c  entries in initial stiffness to rotated versions, then perform r(T) (k r)
c  multiplications.  I need to determine how overwriting occurs in BLAS to
c  determine whether it will then be necessary to copy the new rotated entries
c  back into the original stiffness before reusing the array.
c********
c
      do i=1,nen
        k=ien(i)
        if((skew(1,k).ne.zero).or.(skew(nskdim,k).ne.zero)) then
          call formrt(skew(1,k),rot)
	  call dcopy(nee*nee,s,ione,stemp,ione)
          ll=ndof*(i-1)
c***********  this part still needs to be fixed for BLAS.
          do j=1,nee
            do k=1,ndof
              sum=zero
              do l=1,ndof
                sum=sum+stemp(j,ll+l)*rot(l,k)
              end do
              s(j,ll+k)=sum
            end do
          end do
          call transp(rot,ithree)
	  call dcopy(nee*nee,s,ione,stemp,ione)
          do j=1,ndof
            do k=1,nee
              sum=zero
              do l=1,ndof
                sum=sum+rot(j,l)*stemp(ll+l,k)
              end do
              s(ll+j,k)=sum
            end do
          end do
c**********
        end if
      end do
      return
      end
c
c version
c $Id: rstiff.f,v 1.3 2005/06/08 21:48:12 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
