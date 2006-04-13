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
      subroutine lnklst(neq,lm,lmx,numelv,nen,nee,indx,link,nbrs,iwork,
     & nsizea,nnz,numsn,ierr,errstrng)
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
      integer neq,numelv,nen,nee,iwork,nsizea,nnz,numsn,ierr
      integer lm(nee,numelv),lmx(nee,numelv)
      integer indx(neq),link(iwork),nbrs(iwork)
      character errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer inext,ielg,neeq,i,irow,loc,j,icol
      integer loctmp,ival,locsav
c
cdebug      write(6,*) "Hello from lnklst_f!"
c
      call ifill(indx,izero,neq)
      call ifill(link,izero,iwork)
      call ifill(nbrs,izero,iwork)
c
      neeq=nee
      if(numsn.ne.izero) neeq=2*nee
      inext=ione
      do ielg=1,numelv
c
c      check that available storage is not exceeded
c
        if((inext+nee*nee).gt.iwork) then
          ierr=300
          write(errstrng,700) inext+nee*nee-iwork
          return
        end if
c
        do i=1,neeq
          if(i.le.nee) then
            irow=lm(i,ielg)
          else
            irow=abs(lmx(i-nee,ielg))
          end if
          if(irow.ne.izero) then
            loc=indx(irow)
            do j=1,neeq
              if(j.le.nee) then
                icol=lm(j,ielg)
              else
                icol=abs(lmx(j-nee,ielg))
              end if
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
c $Id: lnklst.f,v 1.6 2005/03/25 23:04:26 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
