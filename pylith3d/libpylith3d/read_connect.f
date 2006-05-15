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
      subroutine read_connect(ien,mat,nen,numelv,numnp,nvfamilies,
     & kr,ifile,ierr,errstrng)
c
c      this subroutine reads element types and connectivities, material
c      types, and infinite element info.
c
c      Error codes:
c          0:  No error
c          1:  Error opening input file
c          3:  Read error
c        106:  Illegal element type
c        107:  Illegal material type
c        108:  Undefined node used for connectivity
c
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nen,numelv,numnp,nvfamilies,kr,ierr
      integer ien(nen,numelv),mat(numelv)
      character ifile*(*),errstrng*(*)
c
c...  local variables
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n,ietypev,infin
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from read_connect_f!"
c
      call ifill(ien,izero,nen*numelv)
      ierr=izero
c
c...  read connectivity, material number, and infinite element info
c
      open(kr,file=ifile,status="old",err=20)
      call pskip(kr)
      do i=1,numelv
        read(kr,*,end=30,err=30) n,ietypev,mat(n),infin,
     &   (ien(j,n),j=1,nen)
c
c...  check for illegal element type
c
        if(ietypev.le.izero.or.ietypev.gt.netypesi) then
          ierr=106
          errstrng="read_connect"
          return
        end if
clater        call infcmp(ietypev,infiel(3,i),inf)
c
c...  check for illegal material type or material model
c
        if(mat(n).le.izero.or.mat(n).gt.nvfamilies) then
          ierr=107
          errstrng="read_connect"
          return
        end if
      end do
      close(kr)
c
c......test for zero or out-of-bound entries in ien array
c
      do n=1,numelv
        do j=1,nen
          if(ien(j,n).le.izero.or.ien(j,n).gt.numnp) then
            ierr=108
            errstrng="read_connect"
            return
          end if
        end do
      end do
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        errstrng="read_connect"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_connect"
        close(kr)
        return
c
      end
c
c version
c $Id: read_connect.f,v 1.11 2005/04/13 00:37:51 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
