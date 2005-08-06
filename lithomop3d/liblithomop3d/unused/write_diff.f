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
      subroutine write_diff(diforc,nslip,idhist,numslp,numdif,numnp,
     & kw,idout,ofile,ierr,errstrng)
c
c...  prints differential forces applied to slippery nodes
c
c     Error codes:
c         0:  No error
c         2:  Error opening output file (if numdif.ne.zero)
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
      integer numslp,numdif,numnp,kw,idout,ierr
      integer nslip(nsdim,numslp),idhist(numnp)
      double precision diforc(ndof,numnp)
      character ofile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer nlines,npage,i,j
      logical nonzed
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
c...  open input file
c
      ierr=izero
      if(numslp.eq.izero.or.numdif.eq.izero.or.idout.eq.izero) return
c
c...  loop over differential forces and output results, if desired
c
      nlines=izero
      npage=50
      open(kw,file=ofile,err=40,status="old",access="append")
      do i=1,numnp
        nonzed=.false.
        do j=1,ndof
          if(diforc(j,i).ne.zero) nonzed=.true.
        end do
        if(nonzed) then
          nlines=nlines+1
          if(nlines.eq.ione.or.mod(nlines,npage).eq.izero) then
            write(kw,6000) (labeld(j),j=1,ndof)
            write(kw,*) ' '
          end if
          write(kw,7000,err=50) i,idhist(i),(diforc(j,i),j=1,ndof)
        end if
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
        errstrng="write_diff"
        close(kw)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="write_diff"
        close(kw)
        return
c
 6000 format(//' differential forces on slippery nodes'//
     & 1x,'  node #   hfac  ',7x,6(a4,11x))
 7000 format(1x,i7,3x,i7,3x,6(3x,1pe12.5))
      end
c
c version
c $Id: write_diff.f,v 1.1 2005/08/05 19:58:07 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
