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
      subroutine printd(d,deld,deltp,kout,ndof,numnp,nout,iflag,idout,
     & idsk,kto,kw,kp)
c
c...program to print displacements
c
c      options:
c
c       iflag=1, displacements are to be output
c       iflag=2, slippery node displacements are to be output using
c                entries in kout to determine which nodes to use.
c
c*  Note:  for full output, call this routine with nout = numnp.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ndof,numnp,nout,iflag,idout,idsk,kto,kw,kp
      integer kout(*)
      double precision d(ndof,numnp),deld(ndof,numnp),deltp
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  defined constants
c
      include "rconsts.inc"
c
      character head(2)*53
      data head/' d i s p l a c e m e n t s                           ',
     &          ' d i f f e r e n t i a l    d i s p l a c e m e n t s'/
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer npage,n,i,m,j
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
cdebug      write(6,*) "Hello from printd_f!"
c
      if(nout.eq.0) return
      npage=50
      if(idout.gt.1.or.idsk.eq.0.or.iflag.eq.2) then
        do n=1,nout
          if((n.eq.1.or.mod(n,npage).eq.0).and.idout.gt.1) then
            write(kw,1000) head(iflag),(labeld(i),i=1,ndof)
            write(kw,*) ' '
          end if
          m=n
          if(iflag.eq.2) m=kout(n)
          if(idout.gt.1) write(kw,2000)  m,(d(i,m),i=1,ndof)
          if(deltp.eq.zero) then
            if(idsk.eq.0) write(kp,3000) m,(d(i,m),i=1,ndof)
            if(idsk.eq.1.and.iflag.eq.2) write(kp) m,(d(i,m),i=1,ndof)
          else
            if(idsk.eq.0) write(kp,3000) m,(d(i,m),i=1,ndof),
     &       (deld(j,m)/deltp,j=1,ndof)
            if(idsk.eq.1.and.iflag.eq.2) write(kp) m,
     &       (d(i,m),i=1,ndof),(deld(j,m)/deltp,j=1,ndof)
          end if
        end do
      end if
      if(idsk.eq.1.and.iflag.eq.1) then
        write(kp) d
        if(deltp.ne.zero) then
          write(kp) deltp
          write(kp) deld
        end if
      end if
      return
1000  format(///1x,a53///4x,'node number',15x,6(a4,14x))
2000  format(6x,i7,10x,6(1pe18.5))
3000  format(i7,12e15.7)
      end
c
c version
c $Id: printd.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
