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
      subroutine scan_coords(nsd,numnp,kr,ierr,coord_units,cfile)
c
c...  subroutine to perform an initial scan of the nodal coordinates
c     file to determine the number of nodal points and the units of
c     the coordinates.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         2:  Units not specified
c         3:  Read error
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nsd,numnp,kr,ierr
      character coord_units*(*),cfile*(*)
c
c...  defined constants
c
      character def(1)*11
      data def/"coord_units"/
c
c...  local variables
c
      integer nget,j,n
      double precision x(3)
      character units(1)*80
      logical units_defined(1)
c
c...  open input file
c
      ierr=0
      numnp=0
      nget=1
      open(kr,file=cfile,status="old",err=20)
c
c...  get coordinate units, returning error 2 if they aren't found.
c
      call get_units(kr,ierr,nget,units_defined,units,def)
      if(ierr.eq.2) return
      coord_units=units(1)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list of coordinates.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) n,(x(j),j=1,nsd)
        numnp=numnp+1
        go to 40
c
c...  normal return
c
 10   continue
        close(kr)
        return
c
c...  error opening file
c
 20   continue
        close(kr)
	ierr=1
        return
c
c...  read error
c
 30   continue
        ierr=3
        close(kr)
        return
c
      end
c
c version
c $Id: scan_coords.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
