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
      subroutine makemsr(ja,indx,link,nbrs,neq,nnz,iwork,nmin,nmax,
     & wavg)
c
c      program to transform linked list into modified sparse row format
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer neq,nnz,iwork,nmin,nmax
      integer ja(nnz),indx(neq),link(iwork),nbrs(iwork)
      double precision wavg
c
c...  defined constants
c
      include "nconsts.inc"
c
c...  intrinsic functions
c
      intrinsic min,max,dble
c
c...  local variables
c
      integer i,loc,ncol
c
      call ifill(ja,izero,nnz)
      ja(1)=neq+2
      nmin=1
      nmax=1
      wavg=neq
      do i=1,neq
        loc=indx(i)
        ncol=0
 20     continue
        ja(ja(i)+ncol)=nbrs(loc)
        ncol=ncol+1
        loc=link(loc)
        if(loc.gt.0) goto 20
        nmin=min(nmin,ncol)
        nmax=max(nmax,ncol)
        wavg=wavg+dble(ncol)
        ja(i+1)=ja(i)+ncol
c
c...    sort entries in each row
c
        call isort(ja(i+1)-ja(i),ja(ja(i)))
      end do
      nmin=nmin+1
      nmax=nmax+1
      if(neq.ne.0) wavg=wavg/dble(neq)
      return
      end
c
c version
c $Id: makemsr.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
