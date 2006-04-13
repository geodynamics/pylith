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
      subroutine id_split(nfault,idftn,numnp,numfn,numflt)
c
c...  subroutine to store the node number of each split node in idftn.
c     This routine needs to be called after nfault has been read
c     (read_split)
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer numnp,numfn,numflt
      integer nfault(3,numfn),idftn(numflt)
c
c...  local variables
c
      integer n,i,nn
c
c... find split nodes and store them in idftn
c
      nn=0
      do n=1,numnp
        do i=1,numfn
          if(n.eq.nfault(2,i)) then
            nn=nn+1
            idftn(nn)=n
            go to 10
          end if
        end do
 10     continue
      end do
c
      return
      end
c
c version
c $Id: id_split.f,v 1.3 2005/04/14 00:56:09 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
