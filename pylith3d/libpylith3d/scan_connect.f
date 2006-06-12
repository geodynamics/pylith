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
      subroutine scan_connect(neni,infmatmod,infmat,ivflist,
     & maxvfamilies,numat,numelv,nvfamilies,ietypev,kr,cfile,ierr,
     & errstrng)
c
c...  subroutine to perform an initial scan of the element connectivity
c     file to determine the total number of elements and the number of
c     elements of each material type.  At present, each material type
c     is assumed to represent an 'element family', and only a single
c     element type is allowed.  In the near future, families will be
c     defined by material model (not type) and element type.
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
      integer maxvfamilies,numat,numelv,nvfamilies,ietypev,kr,ierr
      integer neni(netypesi),infmatmod(6,nmatmodmax),infmat(numat)
      integer ivflist(3,maxvfamilies)
      character cfile*(*),errstrng*(*)
c
c...  local variables
c
      integer j,n,ietype,ietypel,imat,infin
      integer ien(20)
      integer idb
c
      write(6,*) "Hello from scan_connect_f!"
c
c...  open input file
c
      write(6,*) "maxvfamilies:",maxvfamilies
      ierr=izero
      numelv=izero
      ietypev=izero
ctest      nvfamilies=izero
      nvfamilies=maxvfamilies
      call ifill(ivflist,izero,ithree*maxvfamilies)
      open(kr,file=cfile,status="old",err=20)
      do j=1,nvfamilies
        ivflist(2,j)=j
        ivflist(3,j)=infmat(j)
      end do
c
c... scan the file, counting the number of entries for each type of
c    material (element family).
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
      read(kr,*,end=10,err=30) n,ietypel,imat,infin,
     & (ien(j),j=1,neni(ietypel))
      call infcmp(ietypel,ietypev,infin)
      backspace(kr)
      if(ietypel.le.izero.or.ietypel.gt.netypesi) then
        ierr=106
        errstrng="scan_connect"
        return
      end if
 40   continue
        read(kr,*,end=10,err=30) n,ietype,imat,infin,
     &   (ien(j),j=1,neni(ietype))
c
        if(ietype.ne.ietypel) then
          ierr=106
          errstrng="scan_connect"
          return
        end if
c
ctest        if(imat.le.izero.or.imat.gt.numat) then
        if(imat.le.izero) then
          ierr=107
          errstrng="scan_connect"
          return
        end if
c
        ivflist(1,imat)=ivflist(1,imat)+ione
        ivflist(2,imat)=imat
        ivflist(3,imat)=infmat(imat)
        numelv=numelv+ione
        write(6,"(a30)") cfile
cdebug  write(6,*) "n,imat,ivflist:",n,imat,(ivflist(idb,imat),idb=1,3)
c
        go to 40
c
c...  normal return
c
 10   continue
        close(kr)
c
c...  determine the number of element families
c
        do j=1,maxvfamilies
          write(6,*) "j,ivflist(1,j),ivflist(2,j):"
          write(6,*) j,ivflist(1,j),ivflist(2,j)
ctest          if(ivflist(1,j).ne.izero) nvfamilies=nvfamilies+1
        end do
        write(6,*) "nvfamilies:",nvfamilies
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
c $Id: scan_connect.f,v 1.6 2005/04/01 23:24:41 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
