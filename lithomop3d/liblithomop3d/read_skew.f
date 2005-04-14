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
      subroutine read_skew(skew,runits,numrot,numnp,iskopt,kr,sfile,
     & ierr,errstrng)
c
c        program to read and print the array "skew" that specifies
c        local coordinate rotations at each node.  used for applying
c        boundary conditions in directions other than the global
c        coordinate system.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if numrot.ne.zero and
c             iskopt.eq.1)
c         3:  Read error
c         5:  Units not specified
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
      integer numrot,numnp,iskopt,kr,ierr
      double precision skew(nskdim,numnp)
      double precision runits
      character sfile*(*),errstrng*(*)
c
c...  local constants
c
c
c...  intrinsic functions
c
      intrinsic index
c
c...  local variables
c
      integer i,j,n
      character dummy*80
c
      ierr=0
      call fill(skew,zero,nskdim*numnp)
      if(numrot.eq.izero) then
        if(iskopt.eq.itwo) numrot=numnp
        return
      end if
c
c...  open input file and skip over units definitions
c
      open(kr,file=sfile,status="old",err=20)
      call pskip(kr)
      read(kr,"(a80)") dummy
      j=index(dummy,"=")
      if(j.eq.izero) then
        ierr=5
        errstrng="read_skew"
        return
      end if
c
c...  read skew BC
c
      call pskip(kr)
      do i=1,numrot
        read(kr,*,err=30,end=30) n,(skew(j,n),j=1,nskdim)
        do j=1,nskdim
          skew(j,n)=runits*skew(j,n)
        end do
      end do
      close(kr)
c
c...  normal return
c
      return
c
c...  error opening input file
c
20    continue
        ierr=1
        errstrng="read_skew"
        close(kr)
        return
c
c...  error reading input file
c
30    continue
        ierr=3
        errstrng="read_skew"
        close(kr)
        return
c
      end
c
c version
c $Id: read_skew.f,v 1.3 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
