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
      subroutine scan_connect(neni,numeli,netypesi,numelt,nconsize,kr,
     & ierr,cfile)
c
c...  subroutine to perform an initial scan of the element connectivity
c     file to determine the number of each type of element.
c
c     Note that the variable netypesi (and the associated arrays) only
c     defines the number of 'primitive' element types, and does not
c     include the infinite element permutations.  This makes input much
c     simpler, since there are then only 10 element types to define.
c
c     Exceptions should be raised for all errors except code 2, which is
c     not applicable.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         2:  Units not specified (not applicable for this routine)
c         3:  Read error
c         4:  Illegal element type
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer netypesi,numelt,nconsize,kr,ierr
      integer neni(netypesi),numeli(netypesi)
      character cfile*(*)
c
c...  defined constants
c
      include "nconsts.inc"
c
c...  local variables
c
      integer j,n,ietypei,mat,infin
      integer ien(20)
c
c...  open input file
c
      ierr=izero
      numelt=izero
      nconsize=izero
      call ifill(numeli,izero,netypesi)
      open(kr,file=cfile,status="old",err=20)
c
c... scan the file, counting the number of entries for each type of
c    element
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) n,ietypei,mat,infin,
     &   (ien(j),j=1,neni(ietypei))
        if(ietypei.le.izero.or.ietypei.gt.netypesi) then
          ierr=4
          return
        end if
        numeli(ietypei)=numeli(ietypei)+1
        go to 40
c
c...  normal return
c
 10   continue
        do j=1,netypesi
          numelt=numelt+numeli(j)
          nconsize=nconsize+neni(j)*numeli(j)
        end do
        close(kr)
        return
c
c...  error opening file
c
 20   continue
        close(kr)
	ierr=1
        return
c
c...  read error
c
 30   continue
        ierr=3
        close(kr)
        return
c
      end
c
c version
c $Id: scan_connect.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
