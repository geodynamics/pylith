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
      integer inext,iel,indlm,ietype,nee,neeq,i,ii,irow,loc,j,jj,icol
      integer loctmp,ival,locsav
      integer idb,jdb,idbind
c
      write(6,*) "Hello from lnklst_f!"
c
      write(6,*) "neq,nconsz,numelt,iwork,nsizea,nnz,numsn:",neq,nconsz,
     & numelt,iwork,nsizea,nnz,numsn
cdebug      open(15,file="lnklst.info")
cdebug      write(15,*) (lm(idb),idb=1,ndof*nconsz)
cdebug      write(15,*) (lmx(idb),idb=1,ndof*nconsz)
cdebug      do idb=1,numelt
cdebug        write(15,*) (infiel(jdb,idb),jdb=1,6)
cdebug      end do
cdebug      do idb=1,netypes
cdebug        write(15,*) (infetype(jdb,idb),jdb=1,4)
cdebug      end do
cdebug      close(15)
      call ifill(indx,izero,neq)
      call ifill(link,izero,iwork)
      call ifill(nbrs,izero,iwork)
      do idb=1,netypes
        write(6,*) "infetype:",(infetype(jdb,idb),jdb=1,4)
      end do
      write(6,*) "lm:",(lm(idb),idb=1,ndof*nconsz)
      write(6,*) "lmx:",(lm(idb),idb=1,ndof*nconsz)
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
          write(6,*) "Uncaught exception?"
          return
        end if
        indlm=ndof*(infiel(1,iel)-ione)+ione
        ietype=infiel(3,iel)
        nee=infetype(4,ietype)
        neeq=nee
        if(numsn.ne.izero) neeq=itwo*nee
        write(6,*) "iel,indlm,ietype,nee,neeq:",
     &   iel,indlm,ietype,nee,neeq
c
        do i=1,neeq
          ii=indlm+i-1
          if(i.le.nee) irow=lm(ii)
          if(i.gt.nee) irow=abs(lmx(ii-nee))
          if(irow.ne.izero) then
            loc=indx(irow)
            do j=1,neeq
              jj=indlm+j-1
              if(j.le.nee) icol=lm(jj)
              if(j.gt.nee) icol=abs(lmx(jj-nee))
              if(icol.ne.izero.and.icol.ne.irow) then
                write(6,*) "iel,i,j,ii,jj,irow,icol,loc,inext:",
     &           iel,i,j,ii,jj,irow,icol,loc,inext
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
      write(6,*) "After loop in lnklst!"
      call flush(6)
      nsizea=inext-ione
      nnz=nsizea+neq+ione
      write(6,*) "nsizea,nnz,neq:",nsizea,nnz,neq
      call flush(6)
      idbind=max(inext-100,1)
      write(6,*) "end of nbrs:",(nbrs(idb),idb=idbind,inext)
      write(6,*) "end of link:",(link(idb),idb=idbind,inext)
      idbind=max(neq-100,1)
      write(6,*) "end of indx:",(indx(idb),idb=idbind,neq)
      call flush(6)
cdebug      open(15,file="makemsr.info")
cdebug      write(15,*) neq,nnz,iwork
cdebug      write(15,*) (indx(idb),idb=1,neq)
cdebug      write(15,*) (link(idb),idb=1,iwork)
cdebug      write(15,*) (nbrs(idb),idb=1,iwork)
cdebug      close(15)
 700  format("lnklst:  Working storage exceeded by ",i7)
      return
      end
c
c version
c $Id: lnklst-orig.f,v 1.1 2004/08/12 14:40:20 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
