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
      subroutine read_slip(nslip,numslp,numsn,numnp,iskopt,kr,slfile,
     & ierr,errstrng)
c
c...  reads data on free slip interfaces
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if numslp.ne.zero)
c         3:  Read error
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
      integer numslp,numsn,numnp,iskopt,kr,ierr
      integer nslip(nsdim,numslp)
      character slfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n
c
c...  included variable definitions
c
c
c...  open input file
c
      ierr=0
      numsn=0
      if(numslp.eq.izero) return
      call ifill(nslip,izero,nsdim*numslp)
      open(kr,file=slfile,status="old",err=20)
      call pskip(kr)
      do i=1,numslp
        read(kr,*,end=30,err=30) (nslip(j,i),j=1,nsdim)
        if(iskopt.ne.ione) then
          if(nslip(3,i).ne.izero) then
            nslip(4,i)=nslip(3,i)
            nslip(5,i)=izero
          else if(nslip(4,i).ne.izero) then
            nslip(3,i)=nslip(4,i)
            nslip(5,i)=izero
          else if(nslip(5,i).ne.izero) then
            nslip(3,i)=nslip(5,i)
            nslip(4,i)=nslip(5,i)
            nslip(5,i)=izero
          end if
        end if
      end do
c
c...  determine number of slippery nodes
c
      do n=1,numnp
        do j=1,numslp
          if(n.eq.nslip(2,j)) then
            numsn=numsn+ione
            go to 40
          end if
        end do
 40     continue
      end do
      close(kr)
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        errstrng="read_slip"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_slip"
        close(kr)
        return
c
      end
c
c version
c $Id: read_slip.f,v 1.6 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
