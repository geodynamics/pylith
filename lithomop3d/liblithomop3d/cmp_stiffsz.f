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
      subroutine cmp_stiffsz(neq,lm,lmx,infiel,nconsz,numelt,infetype,
     & iwork,numsn,ierr,errstrng)
c
c      program to compute an upper bound on the number of nonzero
c      entries in the stiffness matrix.  Duplicate contributions from
c      different elements are not considered.
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
      integer neq,nconsz,numelt,iwork,numsn,ierr
      integer lm(ndof*nconsz),lmx(ndof*nconsz),infiel(7,numelt)
      integer infetype(4,netypes)
      character errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer iel,indlm,ietype,nee,neeq,i,ii,irow,j,jj,icol
c
cdebug      write(6,*) "Hello from cmp_stiffsz_f!"
c
      iwork=izero
      do iel=1,numelt
        indlm=ndof*(infiel(1,iel)-ione)+ione
        ietype=infiel(3,iel)
        nee=infetype(4,ietype)
        neeq=nee
        if(numsn.ne.izero) neeq=itwo*nee
c
        do i=1,neeq
          ii=indlm+i-1
          irow=lm(ii)
          if(i.gt.nee) irow=abs(lmx(ii-nee))
          if(irow.ne.izero) then
            do j=1,neeq
              jj=indlm+j-1
              icol=lm(jj)
              if(j.gt.nee) icol=abs(lmx(jj-nee))
              if(icol.ne.izero.and.icol.ne.irow) iwork=iwork+ione
            end do
          end if
        end do
      end do
      iwork=iwork+neemax*neemax
      return
      end
c
c version
c $Id: cmp_stiffsz.f,v 1.3 2005/02/23 23:52:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
