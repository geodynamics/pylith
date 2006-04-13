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
      subroutine read_hist(histry,times,nhist,lastep,kr,hfile,
     & ierr,errstrng)
c
c       reads load history definitions, constructs load histories
c       and echos load histories to the output file.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if nhist is greater than
c             zero).
c         3:  Read error
c       111:  Times given are out of order
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nhist,lastep,kr,ierr
      double precision histry(nhist,lastep+1),times(lastep+1)
      character hfile*(*),errstrng*(*)
c
c...  local constants
c
c
c...  intrinsic functions
c
      intrinsic abs,int
c
c... local variables
c
      double precision hloadp,time,hload,diff,diffc,dh,dt,dhdt,defval
      integer i,npoints,j,k,kkp,kk
c
      ierr=izero
      if(nhist.eq.izero) return
c
c...  read load histories
c
      open(kr,file=hfile,status="old",err=20)
      call pskip(kr)
      do i=1,nhist
        read(kr,*,end=30,err=30) npoints,defval
        do j=1,lastep+1
          histry(i,j)=defval
        end do
        hloadp=defval
        kkp=1
        do j=1,npoints
          read(kr,*,end=30,err=30) time,hload
c
c...  find time that most closely matches the given time
c     Note that it is possible to specify a time greater than the last
c     time for which computations are performed.
c
          kk=1
          diff=big
          do k=1,lastep+1
            diffc=abs(time-times(k))
            if(diffc.lt.diff) then
              diff=diffc
              kk=k
            end if
          end do
c
c...  assign loads, using linear interpolation if required
c
          histry(i,kk)=hload
          if(j.gt.ione) then
            dt=times(kk)-times(kkp)
            dh=hload-hloadp
            dhdt=dh/dt
            if(dt.le.zero) then
              ierr=111
              errstrng="read_hist"
              return
            end if
            do k=kkp+1,kk-1
              histry(i,k)=hload+(times(k)-times(kkp))*dhdt
            end do
          end if
          hloadp=hload
          kkp=kk
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
 20   continue
        ierr=1
        errstrng="read_hist"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_hist"
        close(kr)
        return
c
      end
c
c version
c $Id: read_hist.f,v 1.3 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
