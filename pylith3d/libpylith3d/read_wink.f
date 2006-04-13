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
      subroutine read_wink(winkdef,wscal,iwinkdef,iwinkid,nwink,
     & nwinke,kr,wfile,ierr,errstrng)
c
c....program for reading and data on winkler restoring forces
c
c          winkdef(ndof,nwinke) = values of winkler restoring spring
c                             constant, force(i,j)=-wink(i,j)*disp(i,j)
c
c          iwinkdef(ndof,nwinke) = application mode:
c                               iwink = 0, no winkler forces,
c                               iwink = 1, applied throuthout computation
c                               iwink = -n, uses load history factor n
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if nwink.ne.zero)
c         3:  Read error
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
      integer nwink,nwinke,kr,ierr
      integer iwinkdef(ndof,nwinke),iwinkid(nwinke)
      double precision winkdef(ndof,nwinke),wscal(3)
      character wfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j
c
      ierr=izero
      call ifill(iwinkid,izero,nwinke)
      call ifill(iwinkdef,izero,ndof*nwinke)
      call fill(winkdef,zero,ndof*nwinke)
c
c...  open input file
c
      if(nwink.eq.izero) return
      open(kr,file=wfile,status="old",err=20)
c
c.......read winkler force data and output results, if desired
c
      call pskip(kr)
      do i=1,nwinke
        read(kr,*,err=30,end=30) iwinkid(i),(iwinkdef(j,i),j=1,ndof),
     &   (winkdef(j,i),j=1,ndof)
        do j=1,ndof
          winkdef(j,i)=wscal(j)*winkdef(j,i)
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
        errstrng="read_wink"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_wink"
        close(kr)
        return
c
      end
c
c version
c $Id: read_wink.f,v 1.5 2005/04/16 00:45:25 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
