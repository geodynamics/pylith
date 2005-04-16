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
      subroutine cmp_stiffsz(neq,lm,lmx,numelv,iwork,numsn,nen,
     & ierr,errstrng)
c
c      program to compute an upper bound on the number of nonzero
c      entries in the stiffness matrix.  Duplicate contributions from
c      different elements are not considered, so this is quite an
c      overestimate.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer neq,numelv,iwork,numsn,nen,ierr
      integer lm(ndof*nen,numelv),lmx(ndof*nen,numelv)
      character errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer iel,nee,neeq,i,irow,j,icol
c
cdebug      write(6,*) "Hello from cmp_stiffsz_f!"
c
      iwork=izero
      nee=ndof*nen
      neeq=nee
      if(numsn.ne.izero) neeq=itwo*nee
      do iel=1,numelv
c
        do i=1,neeq
          if(i.le.nee) then
            irow=lm(i,iel)
          else
            irow=abs(lmx(i-nee,iel))
          end if
          if(irow.ne.izero) then
            do j=1,neeq
              if(j.le.nee) then
                icol=lm(j,iel)
              else
                icol=abs(lmx(j-nee,iel))
              end if
              if(icol.ne.izero.and.icol.ne.irow) iwork=iwork+ione
            end do
          end if
        end do
      end do
      iwork=iwork+neeq*neeq
      return
      end
c
c version
c $Id: cmp_stiffsz.f,v 1.5 2005/04/16 00:39:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
