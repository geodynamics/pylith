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
      subroutine read_winkx(winkx,wxscal,iwinkx,idx,numnp,ndof,nwinkx,
     & nwinkxe,kr,kw,idout,ierr,wxfile,ofile)
c
c...  program for reading and printing data on differential winkler
c     restoring forces
c
c          winkx(nwinkx) = values of winkler restoring spring
c                             constant, force(idx(i))=-winkx(i)*disp(i)
c
c          iwinkx(2,nwinkx) = application mode:
c                             iwinkx = 0, no winkler forces,
c                             iwinkx = 1, applied throughout computation
c                             iwinkx = -n, uses load history factor n
c
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer numnp,ndof,nwinkx,nwinkxe,kr,kw,idout,ierr
      integer iwinkx(2,nwinkx),idx(ndof,numnp)
      double precision winkx(nwinkx),wxscal(3)
      character wxfile*(*),ofile*(*)
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
      integer i,j,n,nwxtot,nnz,nlines,npage,il
      integer iwxtmp(3)
      double precision wxtmp(3)
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
c...  open output file
c
      if(idout.gt.0) open(kw,file=ofile,status="old",access="append")
c
c...  open input file
c
      ierr=0
      if(nwinkx.eq.0) return
      open(kr,file=wxfile,status="old",err=20)
c
c.......read winkler force data and output results, if desired
c
      call pskip(kr)
      nwxtot=0
      nlines=0
      npage=50
      do i=1,nwinkxe
        read(kr,*,err=30,end=30) n,(iwxtmp(j),j=1,ndof),
     &   (wxtmp(j),j=1,ndof)
        nnz=0
        do j=1,ndof
          wxtmp(j)=wxscal(j)*wxtmp(j)
          if(iwxtmp(j).ne.0) then
            nnz=nnz+1
            nwxtot=nwxtot+1
            iwinkx(1,nwxtot)=iwxtmp(j)
            iwinkx(2,nwxtot)=idx(j,n)
            winkx(nwxtot)=wxtmp(j)
          end if
        end do
        if(idout.gt.0.and.nnz.ne.0) then
          if(mod(nlines,npage).eq.0) then
            write(kw,2000)
            write(kw,3000) (labeld(il),il=1,ndof)
            write(kw,2005)
          end if
          nlines=nlines+1
          write(kw,2020) n,(wxtmp(il),iwxtmp(il),il=1,ndof)
        end if 
      end do
      close(kr)
      if(idout.gt.0) close(kw)
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
 2005 format(/)
 2020 format(1x,i7,2x,3(1pe12.5,2x,i7,3x))
 2000 format(///,'  d i f f e r e n t i a l   w i n k l e r   r e s',
     & ' t o r i n g   f o r c e s',//)
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
c $Id: read_winkx.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
