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
      subroutine read_slip(nslip,numslp,numsn,ndof,numnp,iskopt,nsdim,
     & kr,kw,kp,idout,idsk,ierr,slfile,ofile,pfile)
c
c...  reads and prints data on free slip interfaces
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (no exception should be raised
c             in this case since a slippery node file is optional)
c         2:  Units not specified (not applicable for this routine)
c         3:  Read error
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer numslp,numsn,ndof,numnp,iskopt,nsdim,kr,kw,kp,idout,idsk
      integer ierr
      integer nslip(nsdim,numslp)
      character slfile*(*),ofile*(*),pfile*(*)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
c
c...  intrinsic functions
      intrinsic mod
c
c...  local variables
c
      integer npage,i,j,n
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
c...  open input file
c
      ierr=0
      numsn=0
      call ifill(nslip,izero,nsdim*numslp)
      if(numslp.ne.0) open(kr,file=slfile,status="old",err=20)
c
c...  open output files
c
      if(idout.gt.0) open(kw,file=ofile,status="old",access="append")
      if(idsk.eq.0) open(kp,file=pfile,status="old",access="append")
      if(idsk.eq.1) open(kp,file=pfile,status="old",access="append",
     & form="unformatted")
c
c...  read slippery node entries, and output results if desired
c
      if(idsk.eq.0) write(kp,3000) numslp
      if(idsk.eq.1) write(kp) numslp
      if(numslp.ne.0) then
        call pskip(kr)
        npage=50
        do i=1,numslp
          read(kr,*,end=30,err=30) (nslip(j,i),j=1,nsdim)
          if((i.eq.1.or.mod(i,npage).eq.0).and.idout.gt.0) then
            write(kw,1000) (labeld(j),j=1,ndof)
            write(kw,*) ' '
          end if
          if(iskopt.ne.1.and.ndof.gt.2) then
            if(nslip(3,i).ne.0) then
              nslip(4,i)=nslip(3,i)
              nslip(5,i)=0
            else if(nslip(4,i).ne.0) then
              nslip(3,i)=nslip(4,i)
              nslip(5,i)=0
            else if(nslip(5,i).ne.0) then
              nslip(3,i)=nslip(5,i)
              nslip(4,i)=nslip(5,i)
              nslip(5,i)=0
            end if
          end if
          if(idout.gt.0) write(kw,2000) (nslip(j,i),j=1,ndof+2)
          if(idsk.eq.0) write(kp,3000) (nslip(j,i),j=1,ndof+2)
        end do
      end if
      if(idsk.eq.1.and.numslp.ne.0) write(kp) nslip
c
c...  determine number of slippery nodes
c
      do n=1,numnp
        do j=1,numslp
          if(n.eq.nslip(2,j)) then
            numsn=numsn+1
            go to 40
          end if
        end do
 40     continue
      end do
      if(idsk.eq.0) write(kp,"(i7)") numsn
      if(idsk.eq.1) write(kp) numsn
      close(kr)
      if(idout.gt.1) close(kw)
      close(kp)
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=1
        close(kr)
        return
c
 1000 format(//
     & 'i n t e r n a l   f r e e    s l i p  i n t e r f a c e   d a t
     &a'//
     & '  elem',4x,' node',12x,'dof weight'/
     & '   num',4x,'  num ',6(5x,a4))
 2000 format(i7,4x,i7,1x,6(5x,i4))
 3000 format(16i7)
      end
c
c version
c $Id: read_slip.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
