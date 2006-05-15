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
      subroutine read_coords(x,cscale,numnp,kr,cfile,ierr,errstrng)
c
c...  subroutine to read in nodal coordinates and store them in the x
c     array.  Note that the current implementation would allow a 3D
c     input file to be used for a 2D problem, but not vice-versa.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         3:  Read error
c         5:  Units not specified
c       105:  Wrong number of nodes
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer numnp,kr,ierr
      double precision cscale
      double precision x(nsd,numnp)
      character cfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n
      character dummy*80
c
c...  included variable definitions
c
c
c...  open input file and skip past units specification.
c     Note the the current method of checking for units specification
c     is pretty sloppy, assuming that they were already specified during
c     the scan phase.
c
      ierr=izero
      open(kr,file=cfile,status="old",err=20)
      call pskip(kr)
      read(kr,"(a80)") dummy
      i=index(dummy,"=")
      if(i.eq.izero) then
        ierr=ifive
        errstrng="read_coords"
        return
      end if
      call pskip(kr)
c
c...  read the coordinate entries
c
      do i=1,numnp
        call pskip(kr)
        read(kr,*,end=40,err=30) n,(x(j,n),j=1,nsd)
        do j=1,nsd
          x(j,n)=cscale*x(j,n)
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
        errstrng="read_coords"
        close(kr)
        return
c
c...  read error
c
 30   continue
        ierr=3
        errstrng="read_coords"
        close(kr)
        return
c
c...  too few nodes
c
 40   continue
        ierr=105
        errstrng="read_coords"
        close(kr)
        return
c
      end
c
c version
c $Id: read_coords.f,v 1.4 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
