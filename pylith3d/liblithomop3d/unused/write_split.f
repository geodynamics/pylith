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
      subroutine write_split(fault,nfault,numfn,kw,kp,
     & idout,idsk,ofile,pfile,ierr,errstrng)
c
c...  prints data on split nodes
c
c      fault(ndof,numfn)  = amount of splitting for each dof
c      nfault(3,numfn)    = element, node, and history for each entry
c      numfn              = number of split node entries
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
      integer numfn,kw,kp,idout,idsk,ierr
      integer nfault(3,numfn)
      double precision fault(ndof,numfn)
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer i,j,npage
c
c...  open output files
c
      if(idout.gt.izero) open(kw,file=ofile,err=50,status="old",
     & access="append")
      if(idsk.eq.ione) open(kp,file=pfile,err=50,status="old",
     & access="append")
      if(idsk.eq.itwo) open(kp,file=pfile,err=50,status="old",
     & access="append",form="unformatted")
c
c...  output results if requested
c
      if(idsk.eq.ione) write(kp,"(i7)") numfn
      if(idsk.eq.itwo) write(kp) numfn
      if(idout.gt.izero.or.idsk.eq.ione) then
        npage=50
        do i=1,numfn
          if((i.eq.ione.or.mod(i,npage).eq.izero).and.idout.gt.izero)
     &     write(kw,3000,err=60)
          if(idout.gt.izero) write(kw,4000,err=60) (nfault(j,i),j=1,3),
     &     (fault(j,i),j=1,ndof)
          if(idsk.eq.ione) write(kp,5000,err=60) (nfault(j,i),j=1,2)
        end do
      end if
      if(idsk.eq.itwo.and.numfn.ne.izero) write(kp,err=60) nfault
      if(idsk.eq.itwo.and.numfn.ne.izero) write(kp,err=60) fault
      if(idout.gt.izero) close(kw)
      if(idsk.gt.izero) close(kp)
c
c...  normal return
c
      return
c
c...  error opening output file
c
 50   continue
        ierr=2
        errstrng="write_split"
        if(idsk.gt.izero) close(kp)
        if(idout.gt.izero) close(kw)
        return
c
c...  error writing to output file
c
 60   continue
        ierr=4
        errstrng="write_split"
        if(idsk.gt.izero) close(kp)
        if(idout.gt.izero) close(kw)
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
c $Id: write_split.f,v 1.1 2005/08/05 19:58:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
