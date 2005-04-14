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
      subroutine read_winkx(winkx,wxscal,iwinkx,idx,numnp,nwinkx,
     & nwinkxe,kr,wxfile,ierr,errstrng)
c
c...  program for reading and printing data on differential winkler
c     restoring forces
c
c          winkx(nwinkx) = values of winkler restoring spring
c                             constant, force(idx(i))=-winkx(i)*disp(i)
c
c          iwtmpx(2,nwinkx) = application mode:
c                             iwinkx = 0, no winkler forces,
c                             iwinkx = 1, applied throughout computation
c                             iwinkx = -n, uses load history factor n
c
c        After input, the iwtmp array is transformed into the iwink
c        array, which has the following components:
c          iwinkx(1,i)       = nonzero application mode from iwtmp
c          iwinkx(2,i)       = global equation number to which force is
c                             applied

c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if nwinkx.ne.zero)
c         3:  Read error
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
      integer numnp,nwinkx,nwinkxe,kr,ierr
      integer iwinkx(2,nwinkx),idx(ndof,numnp)
      double precision winkx(nwinkx),wxscal(3)
      character wxfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n,nwxtot,nnz
      integer iwxtmp(3)
      double precision wxtmp(3)
c
c...  included variable definitions
c
c
c...  open output file
c
      ierr=izero
c
c...  open input file
c
      if(nwinkx.eq.izero) return
      open(kr,file=wxfile,status="old",err=20)
c
c.......read winkler force data and output results, if desired
c
      call pskip(kr)
      nwxtot=izero
      do i=1,nwinkxe
        read(kr,*,err=30,end=30) n,(iwxtmp(j),j=1,ndof),
     &   (wxtmp(j),j=1,ndof)
        nnz=izero
        do j=1,ndof
          wxtmp(j)=wxscal(j)*wxtmp(j)
          if(iwxtmp(j).ne.izero) then
            nnz=nnz+1
            nwxtot=nwxtot+1
            iwinkx(1,nwxtot)=iwxtmp(j)
            iwinkx(2,nwxtot)=idx(j,n)
            winkx(nwxtot)=wxtmp(j)
          end if
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
        errstrng="read_winkx"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_winkx"
        close(kr)
        return
c
      end
c
c version
c $Id: read_winkx.f,v 1.3 2005/04/14 00:59:45 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
