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
      subroutine local(id,numnp,ien,lm,infiel,nconsz,numelt,infetype)
c
c.... subroutine to localize id array
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
c
c...  subroutine arguments
c
      integer numnp,nconsz,numelt
      integer id(ndof,numnp),ien(nconsz),lm(ndof,nconsz)
      integer infiel(7,numelt),infetype(4,netypes)
c
c...  local variables
c
      integer iel,j,i,nn,indien,ietype,nen
c
      do iel=1,numelt
        indien=infiel(1,iel)
        ietype=infiel(3,iel)
        nen=infetype(2,ietype)
        do j=indien,indien+nen-1
          nn=ien(j)
          do i=1,ndof
            lm(i,j)=id(i,nn)
          end do
        end do
      end do
      return
      end
c
c version
c $Id: local.f,v 1.3 2005/02/24 00:03:56 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
