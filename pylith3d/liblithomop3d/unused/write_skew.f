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
      subroutine write_skew(skew,numrot,iskopt,numnp,kw,idout,ofile,
     & ierr,errstrng)
c
c        program to print the array "skew" that specifies
c        local coordinate rotations at each node.  used for applying
c        boundary conditions in directions other than the global
c        coordinate system.
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
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numrot,iskopt,numnp,kw,idout,ierr
      double precision skew(nskdim,numnp)
      character ofile*(*),errstrng*(*)
c
c...  local constants
c
      character label(2)*5
      data label/'alpha',' beta'/
c
c...  intrinsic functions
c
      intrinsic index
c
c...  local variables
c
      integer i,j
      logical nonzed
c
      ierr=0
      if(numrot.eq.izero.or.idout.eq.izero) return
c
c...  write skew BC to ascii output file, if desired
c
      open(kw,file=ofile,err=40,status="old",access="append")
      write(kw,1000,err=50) iskopt,(label(i),i=1,nskdim)
      write(kw,1050,err=50)
      do i=1,numnp
        nonzed=.false.
        do j=1,nskdim
          if(skew(j,i).ne.zero) nonzed=.true.
        end do
        if(nonzed) write(kw,2000,err=40) i,(skew(j,i),j=1,nskdim)
      end do
      close(kw)
c
c...  normal return
c
      return
c
c...  error opening output file
c
40    continue
        ierr=2
        errstrng="write_skew"
        close(kw)
        return
c
c...  error writing to output file
c
50    continue
        ierr=4
        errstrng="write_skew"
        close(kw)
        return
c
 1000 format(1x,///,
     &' l o c a l    c o o r d i n a t e    r o t a t i o n s',///,
     &' Coordinate rotation option (iskopt) . . . . . .',i5,/,
     &'     1 = rotations are assigned',/,
     &'     2 = rotations are computed internally',///,
     &'     alpha = cc rotation (radians) from x axis in xy plane',/,
     &'      beta = cc rotation (radians) from x axis in xz plane',//,
     &'   alpha rotation is applied first, then the beta rotation--',/,
     &'   this rotates the skew coordinates into the global system.',//,
     & 1x,'node #',2(10x,a5))
 1050 format(/)
 2000 format(2x,i7,2(f15.2))
      end
c
c version
c $Id: write_skew.f,v 1.1 2005/08/05 19:58:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
