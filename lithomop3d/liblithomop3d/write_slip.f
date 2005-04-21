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
      subroutine write_slip(nslip,numslp,numsn,kw,kp,idout,idsk,
     & ofile,pfile,ierr,errstrng)
c
c...  prints data on free slip interfaces
c
c     Error codes:
c         0:  No error
c         2:  Error opening output file
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
      integer numslp,numsn,kw,kp,idout,idsk,ierr
      integer nslip(nsdim,numslp)
      character ofile*(*),pfile*(*),errstrng*(*)
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
      integer npage,i,j
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
c...  open input file
c
      ierr=0
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
c...  loop over slippery node entries, and output results if desired
c
      if(idsk.eq.ione) write(kp,3000,err=60) numslp
      if(idsk.eq.itwo) write(kp,err=60) numslp
      if(numslp.ne.izero) then
        npage=50
        do i=1,numslp
          if((i.eq.ione.or.mod(i,npage).eq.izero).and.
     &     idout.gt.izero) then
            write(kw,1000,err=60) (labeld(j),j=1,ndof)
            write(kw,*,err=60) ' '
          end if
          if(idout.gt.izero) write(kw,2000,err=60) (nslip(j,i),
     &     j=1,ndof+2)
          if(idsk.eq.ione) write(kp,3000,err=60) (nslip(j,i),
     &     j=1,ndof+2)
        end do
      end if
      if(idsk.eq.itwo.and.numslp.ne.izero) write(kp) nslip
      if(idsk.eq.ione) write(kp,"(i7)",err=60) numsn
      if(idsk.eq.itwo) write(kp,err=60) numsn
      if(idout.gt.izero) close(kw)
      if(idsk.gt.izero) close(kp)
c
c...  normal return
c
      return
c
c...  error opening output file
c
 50   continue
        ierr=2
        errstrng="write_slip"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
        return
c
c...  error writing to output file
c
 60   continue
        ierr=4
        errstrng="write_slip"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
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
c $Id: write_slip.f,v 1.2 2005/04/21 01:07:25 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
