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
      subroutine read_fuldat(iprint,icontr,icode,ncycle,lastep,kr,kw,kp,
     & idout,idsk,fdfile,ofile,pfile,ierr,errstrng)
c
c........reads data on time steps where full outputs are desired
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if there is a
c             time-dependent solution)
c         2:  Error opening output file
c         3:  Read error
c         4:  Write error
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer icontr,icode,ncycle,lastep,kr,kw,kp,idout,idsk,ierr
      integer iprint(icontr)
      character fdfile*(*),ofile*(*),pfile*(*),errstrng*(*)
c
c...  local variables
c
      integer i1,nline,nrem,i,j
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
c...  output results, if desired
c
      if(idout.gt.izero) open(kw,file=ofile,err=40,status="old",
     & access="append")
      if(idsk.eq.izero) open(kp,file=pfile,err=40,status="old",
     & access="append")
      if(idsk.eq.1) open(kp,file=pfile,err=40,status="old",
     & access="append",form="unformatted")
      if(icode.ne.itwo.and.idout.gt.izero) write(kw,1000,err=50)  ncycle
      i1=0
      if((icode.eq.itwo).or.(icontr.eq.izero)) then
        if(idsk.eq.izero) write(kp,4000,err=50) i1,i1,i1
        if(idsk.eq.ione) write(kp,err=50) i1,i1,i1
        close(kp)
        return
      else
        if(idsk.eq.izero) write(kp,4000,err=50) icontr,ncycle,lastep
        if(idsk.eq.ione) write(kp,err=50) icontr,ncycle,lastep
        if(idout.gt.izero) write(kw,2000,err=50) icontr
        if(idsk.eq.ione) write(kp,err=50) (iprint(i),i=1,icontr)
        nline=icontr/7
        nrem=icontr-nline*7
        do i=1,nline
          if(idout.gt.izero) write(kw,3000,err=50) (iprint(7*(i-1)+j),
     &     j=1,7)
          if(idsk.eq.izero) write(kp,4000,err=50) (iprint(7*(i-1)+j),
     &     j=1,7)
        end do
        if(nrem.ne.izero) then
          if(idout.gt.izero) write(kw,3000,err=50) (iprint(7*nline+j),
     &     j=1,nrem)
          if(idsk.eq.izero) write(kp,4000,err=50) (iprint(7*nline+j),
     &     j=1,nrem)
        end if
      end if
      if(idout.gt.izero) close(kw)
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
c...  error opening output file
c
 40   continue
        ierr=2
        errstrng="read_fuldat"
        close(kw)
        close(kp)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="read_fuldat"
        close(kw)
        close(kp)
        return
c
1000  format(//5x,
     &' number of cycles of time sequence . . . . . (ncycle) =',i5/)
2000  format(5x,
     &' the number of full outputs is. . . . . . . .(icontr) =',i5//5x,
     &' full outputs occur at the following timesteps:'//)
3000  format(5x,7(i5,5x))
4000  format(16i5)
c
      end
c
c version
c $Id: read_fuldat.f,v 1.2 2004/07/12 14:35:37 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
