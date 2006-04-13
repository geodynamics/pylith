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
      subroutine nfind(x,xtmp,idslp,ipslp,itmp,itmp1,itmp2,nslip,
     & numslp,numsn,numnp)
c
c       subroutine to find the nsd+1 closest points for each point in
c       idslp.  Also breaks the sets of slippery nodes into separate
c       faults.
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
      integer numslp,numsn,numnp
      integer idslp(numsn),ipslp(npdim,numsn),itmp(numsn),itmp1(numsn)
      integer itmp2(numsn),nslip(nsdim,numslp)
      double precision x(nsd,numnp),xtmp(numsn)
c
c...  intrinsic functions
c
      intrinsic sqrt
c
c...  local variables
c
      integer iflip,isp,ii,n,node,i,nssurf,isn,j,iadd,nn
      double precision xdif
      logical*4 prev
c
c...  loop over number of slippery node entries to find slippery nodes
c     on each fault
c
      iflip=izero
      isp=izero
      ii=izero
      call ifill(itmp1,izero,numsn)
      call ifill(itmp2,izero,numsn)
      do n=1,numslp
        node=nslip(2,n)
        do i=1,ndof
          if(nslip(2+i,n).ne.izero) isn=nslip(2+i,n)
        end do
        if(isn.eq.-isp) iflip=iflip+1
        nssurf=1+iflip/2
        isp=isn
        do i=1,numsn
          if(idslp(i).eq.node) then
            prev=.false.
            do j=1,ii
              if(itmp1(j).eq.node) prev=.true.
            end do
            if(.not.prev) then
              ii=ii+1
              itmp1(ii)=idslp(i)
              itmp2(nssurf)=itmp2(nssurf)+1
            end if
            go to 10
          end if
        end do
10      continue
      end do
      call ifill(ipslp,izero,npdim*numsn)
      iadd=izero
      do nn=1,nssurf
        do n=1,itmp2(nn)
          node=itmp1(n+iadd)
          do i=1,itmp2(nn)
            itmp(i)=i
            xtmp(i)=zero
            do j=1,nsd
              xdif=x(j,itmp1(i+iadd))-x(j,node)
              xtmp(i)=xtmp(i)+xdif*xdif
            end do
            xtmp(i)=sqrt(xtmp(i))
          end do
          call indexx(itmp2(nn),xtmp,itmp)
          do j=1,numsn
            if(idslp(j).eq.node) then
              do i=1,nsd+1
                if(i.le.itmp2(nn)) ipslp(i,j)=itmp1(itmp(i+1)+iadd)
              end do
              go to 20
            end if
          end do
20        continue
        end do
        iadd=iadd+itmp2(nn)
      end do
      return
      end
c
c version
c $Id: nfind.f,v 1.3 2005/04/16 00:43:31 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
