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
      subroutine read_fuldat(iprint,icontr,icode,ncycle,lastep,kr,
     & fdfile,ierr,errstrng)
c
c........reads data on time steps where full outputs are desired
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if there is a
c             time-dependent solution)
c         3:  Read error
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer icontr,icode,ncycle,lastep,kr,ierr
      integer iprint(icontr)
      character fdfile*(*),errstrng*(*)
c
c...  local variables
c
      integer i
c
      ierr=izero
c
c...  read time steps at which a full output is desired
c
      if(lastep.ne.izero.and.icode.eq.ithree) then
        open(kr,file=fdfile,status="old",err=20)
        call pskip(kr)
        do i=1,icontr
          read(kr,*,err=30,end=30) iprint(i)
        end do
        close(kr)
      end if
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        errstrng="read_fuldat"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_fuldat"
        close(kr)
        return
c
      end
c
c version
c $Id: read_fuldat.f,v 1.4 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
