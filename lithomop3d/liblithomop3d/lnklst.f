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
      subroutine lnklst(lm,lmx,indx,link,nbrs,nee,neq,numel,iwork,
     & nsizea,nnz)
c
c      program to create a linked list of nonzero row and column entries
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nee,neq,numel,iwork,nsizea
      integer lm(nee,numel),lmx(nee,numel),indx(neq),link(iwork)
      integer nbrs(iwork)
c
c...  defined constants
c
      include "nconsts.inc"
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer inext,n,i,irow,loc,j,icol,loctmp,ival,locsav,nnz
c
cdebug      write(6,*) "Hello from lnklst_f!"
cdebug      write(6,*) "nee,neq,numel,iwork,nsizea,nnz: ",nee,neq,numel,iwork,
cdebug     & nsizea,nnz
c
      call ifill(indx,izero,neq)
      call ifill(link,izero,iwork)
      call ifill(nbrs,izero,iwork)
c
      inext=1
      do n=1,numel
c
c      check that available storage is not exceeded
c*     Temporary output to stdout that should be replaced by an
c*     exception.
c
        if(inext.gt.iwork) then
          write(6,700) inext-iwork
          stop
        end if
c
        do i=1,2*nee
          irow=lm(i,n)
          if(i.gt.nee) irow=abs(lmx(i-nee,n))
          if(irow.eq.0) goto 20
          loc=indx(irow)
          do j=1,2*nee
            icol=lm(j,n)
            if(j.gt.nee) icol=abs(lmx(j-nee,n))
cdebug            write(6,*) "i,j,irow,icol,loc:",i,j,irow,icol,loc
            if(icol.eq.0) goto 30
            if(icol.eq.irow) goto 30
            if(loc.eq.0) then
              loc=inext
              indx(irow)=inext
              nbrs(inext)=icol
              link(inext)=-irow
              inext=inext+1
            else
              loctmp=loc
  40          continue
              ival=nbrs(loctmp)
              if(icol.eq.ival) goto 30
              locsav=loctmp
              loctmp=link(loctmp)
              if(loctmp.gt.0) goto 40
              link(locsav)=inext
              nbrs(inext)=icol
              link(inext)=-irow
              inext=inext+1
            end if
  30        continue
          end do
  20      continue
        end do
      end do
      nsizea=inext-1
      nnz=nsizea+neq+1
 700  format("Working storage exceeded by ",i7)
      return
      end
c
c version
c $Id: lnklst.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
