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
      subroutine write_timdat(delt,alfa,utol,ftol,etol,times,maxstp,
     & maxit,ntdinit,lgdef,itmax,nintg,lastep,kw,idout,ofile,
     & ierr,errstrng)
c
c...program to write time step data
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
      integer nintg,lastep,kw,idout,ierr
      integer maxstp(nintg),maxit(nintg),ntdinit(nintg),lgdef(nintg)
      integer itmax(nintg)
      double precision delt(nintg),alfa(nintg),utol(nintg),ftol(nintg)
      double precision etol(nintg),times(lastep+1)
      character ofile*(*),errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic index
c
c...  local variables
c
      integer i,j,ii,nstep,nc
c
cdebug      write(6,*) "Hello from write_timdat_f!"
c
      ierr=izero
      if(idout.eq.izero) return
c
c...  echo input to file
c
      open(kw,file=ofile,err=40,status="old",access="append")
      write(kw,1000,err=50)
      do i=1,nintg
          write(kw,2000,err=50) i,maxstp(i),delt(i),alfa(i),maxit(i),
     &   ntdinit(i),lgdef(i),utol(i),ftol(i),etol(i),itmax(i)
      end do
c
c...write out time-step/time relationship for one cycle
c
      nstep=izero
      write(kw,3000,err=50) nstep,1,times(1)
      do i=2,nintg
        ii=i
        write(kw,5000,err=50)
        do j=1,maxstp(i)
          nstep=nstep+ione
          nc=nstep+ione
          write(kw,4000,err=50) nstep,ii,times(nc)
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
 40   continue
        ierr=2
        errstrng="write_timdat"
        close(kw)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="write_timdat"
        close(kw)
        return
c
1000  format(///,' t i m e   s t e p   i n f o r m a t i o n',//,
     & '   Note:  Time step group #1 is the elastic solution',//,
     & 'group   #      step    alfa  maxit ntdinit lgdef utol',
     & '  ftol  etol  itmax',/,
     & '  #   steps    size',/)
2000  format(i3,1x,i5,2x,1pe10.4,2x,0pf4.2,2x,i5,1x,i5,1x,i5,1x,
     & 1pe7.1,1x,1pe7.1,1x,1pe7.1,1x,i5)
3000  format(//,' time-step/time correspondence:',//,
     1 ' time step #   in   group #            time',//,
     2       3x,i5,11x,i5,7x,1pe15.4,5x,'(Elastic)')
4000  format(3x,i5,11x,i5,7x,1pe15.4)
5000  format(' ')
      end
c
c version
c $Id: write_timdat.f,v 1.1 2005/04/16 00:31:53 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
