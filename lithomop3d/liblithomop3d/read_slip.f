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
      subroutine read_slip(nslip,numslp,numsn,numnp,iskopt,
     & kr,kw,kp,idout,idsk,slfile,ofile,pfile,ierr,errstrng)
c
c...  reads and prints data on free slip interfaces
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if numslp.ne.zero)
c         2:  Error opening output file
c         3:  Read error
c         4:  Write error
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
      integer numslp,numsn,numnp,iskopt,kr,kw,kp,idout,idsk,ierr
      integer nslip(nsdim,numslp)
      character slfile*(*),ofile*(*),pfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
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
      if(numslp.ne.izero) open(kr,file=slfile,status="old",err=20)
c
c...  open output files
c
      if(idout.gt.izero) open(kw,file=ofile,err=50,status="old",
     & access="append")
      if(idsk.eq.ione) open(kp,file=pfile,err=50,status="old",
     & access="append")
      if(idsk.eq.itwo) open(kp,file=pfile,err=50,status="old",
     & access="append",form="unformatted")
c
c...  read slippery node entries, and output results if desired
c
      if(idsk.eq.ione) write(kp,3000,err=60) numslp
      if(idsk.eq.itwo) write(kp,err=60) numslp
      if(numslp.ne.izero) then
        call pskip(kr)
        npage=50
        do i=1,numslp
          read(kr,*,end=30,err=30) (nslip(j,i),j=1,nsdim)
          if((i.eq.ione.or.mod(i,npage).eq.izero).and.
     &     idout.gt.izero) then
            write(kw,1000,err=60) (labeld(j),j=1,ndof)
            write(kw,*,err=60) ' '
          end if
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
          if(idout.gt.izero) write(kw,2000,err=60) (nslip(j,i),
     &     j=1,ndof+2)
          if(idsk.eq.ione) write(kp,3000,err=60) (nslip(j,i),
     &     j=1,ndof+2)
        end do
      end if
      if(idsk.eq.itwo.and.numslp.ne.izero) write(kp) nslip
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
      if(idsk.eq.ione) write(kp,"(i7)",err=60) numsn
      if(idsk.eq.itwo) write(kp,err=60) numsn
      close(kr)
      if(idout.gt.ione) close(kw)
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
c...  error opening output file
c
 50   continue
        ierr=2
        errstrng="read_slip"
        close(kr)
        return
c
c...  error writing to output file
c
 60   continue
        ierr=4
        errstrng="read_slip"
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
c $Id: read_slip.f,v 1.4 2004/08/25 01:12:48 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
