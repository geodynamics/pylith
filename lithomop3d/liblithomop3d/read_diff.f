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
      subroutine read_diff(diforc,nslip,idhist,numslp,numdif,ndof,numnp,
     & nsdim,kr,kw,idout,ierr,difile,ofile)
c
c...  reads and prints differential forces applied to slippery nodes
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (no exception should be raised
c             in this case since a differential force file is optional)
c         2:  Units not specified (not yet used for this routine)
c         3:  Read error
c         4:  Differential force applied to non-slippery node
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer numslp,numdif,ndof,numnp,nsdim,kr,kw,idout,ierr
      integer nslip(nsdim,numslp),idhist(numnp)
      double precision diforc(ndof,numnp)
      character difile*(*),ofile*(*)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  data statement definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer nlines,npage,n,i,j
      logical slipry
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
c...  open input file
c
      ierr=0
      call fill(diforc,zero,ndof*numnp)
      call ifill(idhist,izero,numnp)
      if(numslp.eq.0.or.numdif.eq.0) return
      open(kr,file=difile,status="old",err=20)
c
c...  read differential force info and output results, if desired
c
      call pskip(kr)
      nlines=0
      npage=50
      if(idout.gt.0) open(kw,file=ofile,status="old",access="append")
      do i=1,numdif
        read(kr,*,end=30,err=30) n,idhist(n),(diforc(j,n),j=1,ndof)
        slipry=.false.
        if(idout.gt.0) then
          nlines=nlines+1
          if(nlines.eq.1.or.mod(nlines,npage).eq.0) then
            write(kw,6000) (labeld(j),j=1,ndof)
            write(kw,*) ' '
          end if
          write(kw,7000) n,idhist(n),(diforc(j,n),j=1,ndof)
        end if
        do j=1,numslp
          if(n.eq.nslip(2,j)) slipry=.true.
        end do
        if(.not.slipry) then
          ierr=4
          return
        end if
      end do
      close(kr)
      if(idout.gt.0) close(kw)
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        close(kr)
        return
c
 6000 format(//' differential forces on slippery nodes'//
     & 1x,'  node #   hfac  ',7x,6(a4,11x))
 7000 format(1x,i7,3x,i7,3x,6(3x,1pe12.5))
      end
c
c version
c $Id: read_diff.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
