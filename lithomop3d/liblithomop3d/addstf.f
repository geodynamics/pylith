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
      subroutine addstf(alnz,s,lm,lmx,ja,nee,numel,numsn,nnz)
c
c...  subroutine to add local stiffness matrix to global matrix, where
c     the global matrix is stored in modified sparse row (MSR) format.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nee,numel,numsn,nnz
      double precision alnz(nnz),s(nee,nee)
      integer lm(nee*numel),lmx(nee*numel),ja(nnz)
c
c...  intrinsic functions
c
      intrinsic abs,sign,dble
c
c... local variables
c
      integer nloop,i,irow,ii,j,icol,jj,itmp,l
      double precision signi,signj
c
      nloop=nee
      if(numsn.ne.izero) nloop=itwo*nee
      do i=1,nloop
        irow=lm(i)
        signi=one
        ii=i
        if(i.gt.nee) then
          ii=i-nee
          irow=abs(lmx(ii))
          signi=sign(one,dble(lmx(ii)))
        end if
        if(irow.ne.izero) then
          do j=1,nloop
            icol=lm(j)
            signj=one
            jj=j
            if(j.gt.nee) then
              jj=j-nee
              icol=abs(lmx(jj))
              signj=sign(one,dble(lmx(jj)))
            end if
            if(icol.ne.izero) then
              if(irow.eq.icol) then
                alnz(irow)=alnz(irow)+signi*signj*s(ii,jj)
              else
                itmp=izero
                do l=ja(irow),ja(irow+1)-1
                  if(ja(l).eq.icol) then
                    itmp=l
                    go to 10
                  end if
                end do
10              continue
                if(itmp.ne.izero) alnz(itmp)=alnz(itmp)+
     &           signi*signj*s(ii,jj)
              end if
            end if
          end do
        end if
      end do
      return
      end
c
c version
c $Id: addstf.f,v 1.2 2004/06/18 14:50:42 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
