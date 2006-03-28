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
      subroutine write_coords(x,numnp,kw,kp,idout,idsk,ofile,pfile,
     & ierr,errstrng)
c
c...  subroutine to write nodal coordinates to plot and/or ascii files.
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
      integer npage
      parameter(npage=50)
c
c...  subroutine arguments
c
      integer numnp,kw,kp,idout,idsk,ierr
      double precision x(nsd,numnp)
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labelc_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer i,j
c
c...  included variable definitions
c
      include "labelc_def.inc"
c
c...  output coordinates to plot file, if requested
c
      if(idsk.eq.1) then
        call open_append(kp,pfile,"old","",
     &                   ierr,errstrng,"write_coords")
        if(ierr.ne.izero) return
        do i=1,numnp
          write(kp,"(i7,3(1pe16.8))",err=60) i,(x(j,i),j=1,nsd)
        end do
        close(kp)
      else if(idsk.eq.2) then
        call open_append(kp,pfile,"old","unformatted",
     &                   ierr,errstrng,"write_coords")
        if(ierr.ne.izero) return
        write(kp,err=60) x
        close(kp)
      end if
c
c...  output coordinates to human-readable file
c
      if(idout.gt.izero) then
        call open_append(kw,ofile,"old","",
     &                   ierr,errstrng,"write_coords")
        if(ierr.ne.izero) return
        do i=1,numnp
          if(i.eq.1.or.mod(i,npage).eq.0) write(kw,1000)
     &     (labelc(j),j=1,nsd)
          write(kw,"(6x,i7,10x,3(1pe20.8))",err=60) i,(x(j,i),j=1,nsd)
        end do
        close(kw)
      end if
c
c...  normal return
c
      return
c
c...  write error
c
 60   continue
        ierr=4
        errstrng="write_coords"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
        return
c
1000  format(1x,///,' n o d a l   c o o r d i n a t e   d a t a',///,
     & 4x,' node number ',13x,3(a4,18x),//)
c
      end
c
c version
c $Id: write_coords.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
