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
      subroutine read_split(fault,nfault,numfn,numflt,ndof,numnp,kr,kw,
     & kp,idout,idsk,ierr,spfile,ofile,pfile)
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
c         1:  Error opening input file (no exception should be raised
c             in this case since a split node file is optional)
c         2:  Units not specified (not yet used for this routine)
c         3:  Read error
c

      include "implicit.inc"
c
c...  subroutine arguments
c
      integer numfn,numflt,ndof,numnp,kr,kw,kp,idout,idsk,ierr
      integer nfault(3,numfn)
      double precision fault(ndof,numfn)
      character spfile*(*),ofile*(*),pfile*(*)
c
c...  defined constants
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
      integer i,j,n,npage
c
c...  open input file and read entries
c
      ierr=0
      numflt=0
      if(numfn.ne.0) then
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
c...  open output files
c
      if(idout.gt.0) open(kw,file=ofile,status="old",access="append")
      if(idsk.eq.0) open(kp,file=pfile,status="old",access="append")
      if(idsk.eq.1) open(kp,file=pfile,status="old",access="append",
     & form="unformatted")
c
c...  output results if requested
c
      if(idsk.eq.0) write(kp,"(i7)") numfn
      if(idsk.eq.1) write(kp) numfn
      if(idout.gt.0.or.idsk.eq.0) then
        npage=50
        do i=1,numfn
          if((i.eq.1.or.mod(i,npage).eq.0).and.idout.gt.0)
     &     write(kw,3000)
          if(idout.gt.0) write(kw,4000) (nfault(j,i),j=1,3),
     &     (fault(j,i),j=1,ndof)
          if(idsk.eq.0) write(kp,5000) (nfault(j,i),j=1,2)
        end do
      end if
      if(idsk.eq.1.and.numfn.ne.0) write(kp) nfault
      if(idsk.eq.1.and.numfn.ne.0) write(kp) fault
c
c     compute total number of split nodes
c
      if(numfn.ne.0) then
        do n=1,numnp
          do i=1,numfn
            if(n.eq.nfault(2,i)) then
              numflt=numflt+1
              goto 40
            end if
          end do
 40       continue
        end do
      end if
      close(kr)
      if(idout.gt.0) close(kw)
      close(kp)
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
 3000 format(1x///1x,'s p l i t  n o d e   d a t a'//
     &'     ihist =  0, slip in elastic computation'/
     &'           = -1, slip at constant velocity '/
     &'           =  n, uses load history factor n'//
     &' elem',6x,'node',7x,'ihist',7x,'split',12x,'split',12x,'split'/
     &'  num       num',20x,'dof1',13x,'dof2',13x,'dof3'//)
 4000 format(i7,5x,i7,5x,i7,3(5x,1pe12.5))
 5000 format(2i7)
      end
c
c version
c $Id: read_split.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
