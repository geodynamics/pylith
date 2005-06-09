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
