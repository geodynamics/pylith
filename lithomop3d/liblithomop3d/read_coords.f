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
      subroutine read_coords(x,cscale,numnp,kr,kw,kp,idout,idsk,
     & cfile,ofile,pfile,ierr,errstrng)
c
c...  subroutine to read in nodal coordinates and store them in the x
c     array.  Note that the current implementation would allow a 3D
c     input file to be used for a 2D problem, but not vice-versa.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         2:  Error opening output file
c         3:  Read error
c         4:  Write error
c         5:  Units not specified
c       105:  Wrong number of nodes
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
      integer numnp,kr,kw,kp,idout,idsk,ierr
      double precision cscale
      double precision x(nsd,numnp)
      character cfile*(*),ofile*(*),pfile*(*),errstrng*(*)
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
      integer i,j,n
      character dummy*80
c
c...  included variable definitions
c
      include "labelc_def.inc"
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
        read(kr,*,end=40,err=30) n,(x(j,i),j=1,nsd)
        do j=1,nsd
          x(j,i)=cscale*x(j,i)
        end do
      end do
      close(kr)
c
c...  output coordinates to plot file
c
      if(idsk.eq.1) then
        open(kp,file=pfile,err=50,status="old",access="append")
        do i=1,numnp
          write(kp,"(i7,3(1pe16.8))",err=60) i,(x(j,i),j=1,nsd)
        end do
      else if(idsk.eq.2) then
        open(kp,file=pfile,err=50,status="old",access="append",
     &   form="unformatted")
        write(kp,err=60) x
      end if
      close(kp)
c
c...  output coordinates to human-readable file
c
      if(idout.gt.izero) then
        open(kw,file=ofile,err=50,status="old",access="append")
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
c...  error opening output file
c
 50   continue
        ierr=2
        errstrng="read_coords"
        close(kr)
        return
c
c...  write error
c
 60   continue
        ierr=4
        errstrng="read_coords"
        close(kr)
        return
c
1000  format(1x,///,' n o d a l   c o o r d i n a t e   d a t a',///,
     & 4x,' node number ',13x,3(a4,18x),//)
c
      end
c
c version
c $Id: read_coords.f,v 1.3 2004/08/25 01:12:48 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
