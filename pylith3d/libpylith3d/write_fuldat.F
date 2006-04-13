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
      subroutine write_fuldat(iprint,icontr,icode,ncycle,lastep,kw,kp,
     & idout,idsk,ofile,pfile,ierr,errstrng)
c
c...  writes data on time steps where full outputs are desired
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
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer icontr,icode,ncycle,lastep,kw,kp,idout,idsk,ierr
      integer iprint(icontr)
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  local variables
c
      integer i1,nline,nrem,i,j
c
      ierr=izero
c
c...  output results, if desired
c
      if(idout.gt.izero) then
        call open_append(kw,ofile,"old","",ierr,errstrng,"write_fuldat")
        if(ierr.ne.izero) return
      end if
      if(idsk.eq.ione) then
        call open_append(kp,pfile,"old","",ierr,errstrng,"write_fuldat")
        if(ierr.ne.izero) return
      else if(idsk.eq.itwo) then
        call open_append(kp,pfile,"old","unformatted",
     &                   ierr,errstrng,"write_fuldat")
        if(ierr.ne.izero) return
      end if
      if(icode.ne.itwo.and.idout.gt.izero) write(kw,1000,err=50)  ncycle
      i1=0
      if((icode.eq.itwo).or.(icontr.eq.izero)) then
        if(idsk.eq.ione) write(kp,4000,err=50) i1,i1,i1
        if(idsk.eq.itwo) write(kp,err=50) i1,i1,i1
      else
        if(idsk.eq.ione) write(kp,4000,err=50) icontr,ncycle,lastep
        if(idsk.eq.itwo) write(kp,err=50) icontr,ncycle,lastep
        if(idout.gt.izero) write(kw,2000,err=50) icontr
        if(idsk.eq.itwo) write(kp,err=50) (iprint(i),i=1,icontr)
        nline=icontr/7
        nrem=icontr-nline*7
        do i=1,nline
          if(idout.gt.izero) write(kw,3000,err=50) (iprint(7*(i-1)+j),
     &     j=1,7)
          if(idsk.eq.ione) write(kp,4000,err=50) (iprint(7*(i-1)+j),
     &     j=1,7)
        end do
        if(nrem.ne.izero) then
          if(idout.gt.izero) write(kw,3000,err=50) (iprint(7*nline+j),
     &     j=1,nrem)
          if(idsk.eq.ione) write(kp,4000,err=50) (iprint(7*nline+j),
     &     j=1,nrem)
        end if
      end if
      if(idout.gt.izero) close(kw)
      if(idsk.gt.izero) close(kp)
c
c...  normal return
c
      return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="write_fuldat"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
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
c $Id: write_fuldat.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
