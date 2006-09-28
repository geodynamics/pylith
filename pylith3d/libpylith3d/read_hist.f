c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
