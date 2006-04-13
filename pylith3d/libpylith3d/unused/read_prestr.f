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
      subroutine read_prestr(stn,st0,prscal,nstr,ngauss,nppts,numel,
     & nprestr,ipstrs,idout,kr,kw,ierr,prfile,ofile)
c
c.......program to read and print background stress data
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (this should not produce an
c             exception since prestresses are optional)
c         2:  Units not specified (not applicable for this case)
c         3:  Read error
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,ngauss,nppts,numel,nprestr,ipstrs,idout,kr,kw,ierr
      double precision stn(nstr,ngauss,numel),st0(nstr,nppts,numel)
      character prfile*(*),ofile*(*)
c
c...  included dimension and type statements
c
      include "prscal_dim.inc"
      include "labels_dim.inc"
c
c...  defined constants
c
      include "rconsts.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer nlines,interv,i,l,j,n,ll,iline
c
c...  included variable definitions
c
      include "labels_def.inc"
c
c...  return if there are no prestresses or if automatic prestress
c     computation is being performed
c
      ierr=0
      call fill(st0,zero,nstr*nppts*numel)
      if(nprestr.eq.0.or.ipstrs.eq.1) return
c
c...  open input file
c
      open(kr,file=prfile,status="old",err=20)
      if(idout.gt.0) open(kw,file=ofile,status="old",access="append")
c
c...  read and scale prestress entries.
c     Transfer the prestresses into the intial stress values.
c     Output prestresses if requested.
c
      call fill(stn,zero,nstr*ngauss*numel)
      call pskip(kr)
      nlines=0
      interv=50/nppts
      do i=1,nprestr
        do l=1,nppts
          read(kr,*,end=30,err=30) n,(st0(j,l,n),j=1,nstr)
          do j=1,nstr
            st0(j,l,n)=prscal(j)*st0(j,l,n)
          end do
        end do
        do l=1,ngauss
          ll=l
          if(nppts.eq.1) ll=1
          do j=1,nstr
            stn(j,l,n)=st0(j,ll,n)
          end do
        end do
        if(idout.gt.0) then
          if(mod(nlines,interv).eq.0) then
            write(kw,1000) (labels(iline),iline=1,nstr)
            write(kw,'(/)')
          end if
          nlines=nlines+nppts
          do l=1,nppts
            if(idout.gt.0) write(kw,2000) n,(st0(j,l,n),j=1,nstr)
          end do
        end if
      end do
      close(kr)
      if(idout.gt.0) close(kw)
c
c...  normal return
c
      return
c
c...  error opening file
c
 20   continue
        ierr=1
        close(kr)
        return
c
c...  read error
c
 30   continue
        ierr=3
        close(kr)
        return
c
 1000 format(1x///' b a c k g r o u n d   s t r e s s   d a t a'//
     1    3x,'elem #',6(4x,a6,6x))
 2000 format(1x,i7,5x,6(1pe12.5,4x))
      end
c
c version
c $Id: read_prestr.f,v 1.1 2004/07/12 18:17:42 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
