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
      subroutine write_hist(histry,times,nhist,lastep,kw,idout,
     & ofile,ierr,errstrng)
c
c       writes load histories to the output file.
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
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nhist,lastep,kw,idout,ierr
      double precision histry(nhist,lastep+1),times(lastep+1)
      character ofile*(*),errstrng*(*)
c
c...  local constants
c
      character*1 star(30)
      data star/30*'*'/
c
c...  intrinsic functions
c
      intrinsic abs,int
c
c... local variables
c
      double precision fmax,fmin,fac,range
      integer i,j,nstars,l
c
      ierr=izero
      if(idout.eq.izero) return
      call open_append(kw,ofile,"old","",ierr,errstrng,"write_hist")
      if(ierr.ne.izero) return
      write(kw,1000,err=70) nhist
      if(nhist.eq.izero) close(kw)
      if(nhist.eq.izero) return
c
c...  write load histories
c
      do i=1,nhist
        write(kw,2000,err=70) i
c
c     echo load history and construct a miniplot
c
        fmax=-1.d32
        fmin=1.d32
        do j=1,lastep+1
          fac=histry(i,j)
          if(fac.gt.fmax) fmax=fac
          if(fac.lt.fmin) fmin=fac
        end do
        if(fmin.gt.zero) fmin=zero
        range=fmax-fmin
        do j=1,lastep+1
          if(range.ne.0) then
            nstars=int(30.0d0*((histry(i,j)-fmin)/range)+.001d0)
          else
            nstars=30
          end if
          write(kw,3000,err=70) j-1,times(j),histry(i,j),
     &     '|',(star(l),l=1,nstars)
        end do
      end do
      close(kw)
c
c...  normal return
c
      return
c
c...  error opening output file
c
 60   continue
        ierr=2
        errstrng="write_hist"
        close(kw)
        return
c
c...  error writing to output file
c
 70   continue
        ierr=4
        errstrng="write_hist"
        close(kw)
        return
c
1000  format(///' l o a d   h i s t o r y   f a c t o r s'///,
     &          ' number of load histories defined (nhist) ...',i2/)
2000  format(// ' load history factor # ',i2//
     1  '  time step       time            factor',
     2 '               scaled miniplot'//)
3000  format(2x,i5,7x,1pe12.4,4x,1pe12.4,5x,31a1)
      end
c
c version
c $Id: write_hist.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
