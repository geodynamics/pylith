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
      subroutine read_wink(wink,wscal,iwink,id,numnp,nwink,nwinke,
     & kr,kw,idout,wfile,ofile,ierr,errstrng)
c
c....program for reading and printing data on winkler restoring forces
c
c          wink(ndof,numnp) = values of winkler restoring spring
c                             constant, force(i,j)=-wink(i,j)*disp(i,j)
c
c          iwink(ndof,numnp) = application mode:
c                               iwink = 0, no winkler forces,
c                               iwink = 1, applied throuthout computation
c                               iwink = -n, uses load history factor n
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (only if nwink.ne.zero)
c         2:  Error opening output file
c         3:  Read error
c         4:  Write error
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
      integer numnp,nwink,nwinke,kr,kw,idout,ierr
      integer iwink(2,nwink),id(ndof,numnp)
      double precision wink(nwink),wscal(3)
      character wfile*(*),ofile*(*),errstrng*(*)
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
      integer i,j,n,nwtot,nnz,nlines,npage,iline
      integer iwtmp(3)
      double precision wtmp(3)
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
c...  open output file
c
      ierr=izero
      if(idout.gt.izero) open(kw,file=ofile,err=40,status="old",
     & access="append")
c
c...  open input file
c
      if(nwink.eq.izero) return
      open(kr,file=wfile,status="old",err=20)
c
c.......read winkler force data and output results, if desired
c
      call pskip(kr)
      nwtot=izero
      nlines=izero
      npage=50
      do i=1,nwinke
        read(kr,*,err=30,end=30) n,(iwtmp(j),j=1,ndof),
     &   (wtmp(j),j=1,ndof)
        nnz=izero
        do j=1,ndof
          wtmp(j)=wscal(j)*wtmp(j)
          if(iwtmp(j).ne.izero) then
            nnz=nnz+1
            nwtot=nwtot+1
            iwink(1,nwtot)=iwtmp(j)
            iwink(2,nwtot)=id(j,n)
            wink(nwtot)=wtmp(j)
          end if
        end do
        if(idout.gt.izero.and.nnz.ne.izero) then
          if(mod(nlines,npage).eq.izero) then
            write(kw,1000,err=50)
            write(kw,3000,err=50) (labeld(iline),iline=1,ndof)
            write(kw,2005,err=50)
          end if
          nlines=nlines+1
          write(kw,2020,err=50) n,(wtmp(iline),iwtmp(iline),
     &     iline=1,ndof)
        end if 
      end do
      close(kr)
      if(idout.gt.izero) close(kw)
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
c...  error opening output file
c
 40   continue
        ierr=2
        errstrng="read_wink"
        close(kw)
        return
c
c...  error reading input file
c
 50   continue
        ierr=4
        errstrng="read_wink"
        close(kw)
        return
c
 2005 format(/)
 2020 format(1x,i7,2x,3(1pe12.5,2x,i7,3x))
 1000 format(///,'  w i n k l e r   r e s t o r i n g   f o r c e s',//)
c 2000 format(///,'  d i f f e r e n t i a l   w i n k l e r   r e s',
c     & ' t o r i n g   f o r c e s',//)
 3000 format(
     x '      key to application mode codes:',//
     x '          mode = 0, no winkler restoring forces; the',/,
     x '                 numerical value to the left is meaningless.',/,
     x '          mode = 1, forces applied at all times',/,
     x '          mode = -n, uses n^th load history factor',///,
     x 1x,'node # ',a4,' winklr     mode',5(2x,a4,' winklr     mode'),/)
      end
c
c version
c $Id: read_wink.f,v 1.3 2005/04/02 22:11:15 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
