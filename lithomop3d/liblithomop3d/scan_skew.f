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
      subroutine scan_skew(nskdim,numrot,kr,ierr,rotation_units,skfile)
c
c...  subroutine to perform an initial scan of the skew rotation file
c     to determine the number of predefined skew rotations and the
c     units being used for rotations.
c
c     Note that the total number of skew rotations could be greater than
c     the value determined here if automatic rotation computation is being
c     performed for slippery nodes.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (this should not produce an
c             exception since skew BC are optional)
c         2:  Units not specified
c         3:  Read error
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nskdim,numrot,kr,ierr
      character rotation_units*(*),skfile*(*)
c
c...  defined constants
c
      character def(1)*14
      data def/"rotation_units"/
c
c...  local variables
c
      integer nget,j,n
      double precision skew(2)
      character units(1)*80
      logical units_defined(1)
c
c...  open input file
c
      ierr=0
      numrot=0
      nget=1
      open(kr,file=skfile,status="old",err=10)
c
c...  get units, returning error 2 if they aren't found.
c
      call get_units(kr,ierr,nget,units_defined,units,def)
      if(ierr.eq.2) return
      rotation_units=units(1)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) n,(skew(j),j=1,nskdim)
        numrot=numrot+1
        go to 40
c
c...  normal return
c
 10   continue
        close(kr)
        return
c
c...  error opening file -- not used since this file is optional
c
c* 20   continue
c*        ierr=1
c*        close(kr)
c*        return
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
c $Id: scan_skew.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
