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
      subroutine adjid(id,idx,nslip,idslp,numslp,numnp,numsn,neq)
c
c      adjusts id array and creates idx array for additional degrees
c      of freedom associated with free slip interfaces
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer numslp,numnp,numsn,neq
      integer id(ndof,numnp),idx(ndof,numnp),nslip(nsdim,numslp)
      integer idslp(numsn)
c
c...  local variables
c
      integer n,j,i,ii,nn,ns
c
      call ifill(idx,izero,ndof*numnp)
      if(numslp.eq.izero) return
      ns=izero
      do n=1,numnp
        do j=1,numslp
          if(n.eq.nslip(2,j)) then
            ns=ns+1
            idslp(ns)=n
            do i=1,ndof
              if((nslip(2+i,j).ne.izero).and.(id(i,n).ne.izero)) then
                idx(i,n)=id(i,n)
                do ii=i,ndof
                  if(id(ii,n).ne.izero) id(ii,n)=id(ii,n)+1
                end do
                do nn=n+1,numnp
                  do ii=1,ndof
                    if(id(ii,nn).ne.izero) id(ii,nn)=id(ii,nn)+1
                  end do
                end do
                neq=neq+1
              end if
            end do
            goto 100
          end if
        end do
  100   continue
      end do
      return
      end
c
c version
c $Id: adjid.f,v 1.2 2004/06/18 14:55:02 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
