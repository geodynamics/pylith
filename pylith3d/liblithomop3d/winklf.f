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
      subroutine winklf(bwink,bintern,deld,iwink,wink,histry,nwink,
     & nhist,nstep,neq,lastep)
c
c       program to compute winkler restoring forces from the
c       displacements and add them to the internal force vector bintern.
c       Note:  This method differs from the previous method since the
c       Winkler forces are now included as part of the internal force
c       vector.  Thus, the forces are now positive.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nwink,nhist,nstep,neq,lastep
      integer iwink(2,nwink)
      double precision bwink(neq),bintern(neq),deld(neq),wink(nwink)
      double precision histry(nhist,lastep+1)
c
c...  local variables
c
      integer i,mode,k,ihist
c
cdebug      write(6,*) "Hello from winklf_f!"
c
      do i=1,nwink
        mode=iwink(1,i)
        k=iwink(2,i)
        if(mode.eq.ione) then
          bwink(k)=bwink(k)+wink(i)*deld(k)
        else
          ihist=-mode
          bwink(k)=bwink(k)+wink(i)*deld(k)*histry(ihist,nstep+1)
        end if
        bintern(k)=bintern(k)+bwink(k)
      end do
      return
      end
c
c version
c $Id: winklf.f,v 1.4 2005/05/03 18:41:29 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
