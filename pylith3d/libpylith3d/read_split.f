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
      subroutine read_split(fault,nfault,numfn,numflt,numnp,numelv,
     & kr,spfile,ierr,errstrng)
c
c...  reads and prints data on split nodes
c
c      fault(ndof,numfn)  = amount of splitting for each dof
c      nfault(3,numfn)    = element, node, and history for each entry
c      numfn              = number of split node entries
c      numflt             = total number of split nodes
c
c     Once the total number of split nodes has been computed, routine
c     id_split should be called to store the entries in idftn.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if numfn.ne.zero)
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
      integer numfn,numflt,numnp,numelv,kr,ierr
      integer nfault(3,numfn)
      double precision fault(ndof,numfn)
      character spfile*(*),errstrng*(*)
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n
c
c...  open input file and read entries
c
      ierr=izero
      numflt=izero
      if(numfn.ne.izero) then
        call fill(fault,zero,ndof*numfn)
        call ifill(nfault,izero,3*numfn)
        open(kr,file=spfile,status="old",err=20)
        call pskip(kr)
        do i=1,numfn
          read(kr,*,end=30,err=30) (nfault(j,i),j=1,3),
     &     (fault(j,i),j=1,ndof)
        end do
      end if
c
c     compute total number of split nodes
c
      if(numfn.ne.izero) then
        do n=1,numnp
          do i=1,numfn
            if(n.eq.nfault(2,i)) then
              numflt=numflt+ione
              goto 40
            end if
          end do
 40       continue
        end do
      end if
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
        errstrng="read_split"
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_split"
        close(kr)
        return
c
      end
c
c version
c $Id: read_split.f,v 1.5 2005/04/14 00:59:44 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
