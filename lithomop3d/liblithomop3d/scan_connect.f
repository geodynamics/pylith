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
      subroutine scan_connect(neni,infmat,infmatmod,matmod,
     & indprop,numat,numelt,nconsz,kr,cfile,ierr,errstrng)
c
c...  subroutine to perform an initial scan of the element connectivity
c     file to determine the total number of elements and the number of
c     elements of each material type.
c
c     Note that the variable netypesi (and the associated arrays) only
c     defines the number of 'primitive' element types, and does not
c     include the infinite element permutations.  This makes input much
c     simpler, since there are then only 10 element types to define.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         3:  Read error
c       106:  Illegal element type
c       101:  Attempt to use undefined material model
c       107:  Illegal material type
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer numat,numelt,nconsz,kr,ierr
      integer neni(netypesi),infmat(3,numat),infmatmod(5,nmatmodmax)
      integer matmod(numat),indprop(numat)
      character cfile*(*),errstrng*(*)
c
c...  local variables
c
      integer i,j,n,ietypei,imat,infin,matm
      integer ien(20)
c
c...  open input file
c
      ierr=izero
      numelt=izero
      nconsz=izero
      call ifill(infmat,izero,ithree*numat)
      open(kr,file=cfile,status="old",err=20)
c
c...  set up part of the infmat array that has already been determined
c     by python material property input routine.
c
      do i=1,numat
        infmat(1,i)=matmod(i)
        infmat(3,i)=indprop(i)
      end do
c
c... scan the file, counting the number of entries for each type of
c    material.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) n,ietypei,imat,infin,
     &   (ien(j),j=1,neni(ietypei))
c
        matm=infmat(1,imat)
c
        if(ietypei.le.izero.or.ietypei.gt.netypesi) then
          ierr=106
          errstrng="scan_connect"
          return
        end if
c
        if(imat.le.izero.or.imat.gt.numat) then
          ierr=107
          errstrng="scan_connect"
          return
        end if
c
        if(matm.gt.nmatmodmax.or.infmatmod(1,matm).eq.izero) then
          ierr=101
          errstrng="scan_connect"
          return
        end if
c
        numelt=numelt+ione
        nconsz=nconsz+neni(ietypei)
        infmat(2,imat)=infmat(2,imat)+ione
c
        go to 40
c
c...  normal return
c
 10   continue
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
c $Id: scan_connect.f,v 1.2 2004/07/09 17:01:13 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
