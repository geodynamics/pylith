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
      subroutine lnklst(neq,lm,lmx,infiel,nconsz,numelt,infetype,indx,
     & link,nbrs,iwork,nsizea,nnz,numsn,ierr,errstrng)
c
c      program to create a linked list of nonzero row and column entries
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
      integer neq,nconsz,numelt,iwork,nsizea,nnz,numsn,ierr
      integer lm(ndof*nconsz),lmx(ndof*nconsz),infiel(6,numelt)
      integer infetype(4,netypes),indx(neq),link(iwork),nbrs(iwork)
      character errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer inext,iel,indien,ietype,nee,neeq,i,ii,irow,loc,j,jj,icol
      integer loctmp,ival,locsav
c
cdebug      write(6,*) "Hello from lnklst_f!"
c
      call ifill(indx,izero,neq)
      call ifill(link,izero,iwork)
      call ifill(nbrs,izero,iwork)
c
      inext=ione
      do iel=1,numelt
c
c      check that available storage is not exceeded
c*     Temporary output to stdout that should be replaced by an
c*     exception.
c
        if(inext.gt.iwork) then
          ierr=300
          write(errstrng,700) inext-iwork
          return
        end if
        indien=infiel(1,iel)
        ietype=infiel(3,iel)
        nee=infetype(4,ietype)
        neeq=nee
        if(numsn.ne.izero) neeq=itwo*nee
c
        do i=1,neeq
          ii=indien+i-1
          if(i.le.nee) irow=lm(ii)
          if(i.gt.nee) irow=abs(lmx(ii-nee))
          if(irow.ne.izero) then
            loc=indx(irow)
            do j=1,neeq
              jj=indien+j-1
              if(j.le.nee) icol=lm(jj)
              if(j.gt.nee) icol=abs(lmx(jj-nee))
cdebug            write(6,*) "i,j,irow,icol,loc:",i,j,irow,icol,loc
              if(icol.ne.izero.and.icol.ne.irow) then
                if(loc.eq.izero) then
                  loc=inext
                  indx(irow)=inext
                  nbrs(inext)=icol
                  link(inext)=-irow
                  inext=inext+1
                else
                  loctmp=loc
 40               continue
                    ival=nbrs(loctmp)
                    if(icol.eq.ival) goto 30
                    locsav=loctmp
                    loctmp=link(loctmp)
                  if(loctmp.gt.izero) goto 40
                  link(locsav)=inext
                  nbrs(inext)=icol
                  link(inext)=-irow
                  inext=inext+1
                end if
              end if
 30           continue
            end do
          end if
        end do
      end do
      nsizea=inext-ione
      nnz=nsizea+neq+ione
 700  format("lnklst:  Working storage exceeded by ",i7)
      return
      end
c
c version
c $Id: lnklst.f,v 1.3 2004/07/16 16:02:21 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
