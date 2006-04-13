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
      subroutine printl(idx,iwinkx,idslp,histry,numsn,numnp,nstep,
     & nhist,nwinkx,lastep,idsk,kp)
c
c...subroutine to print array that specifies whether each slippery
c   node is locked (1) or unlocked (0) during the current time step.
c   This information is used for postprocessing.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer numsn,numnp,nstep,nhist,nwinkx,lastep,idsk,kp
      integer idx(ndof,numnp),iwinkx(2,nwinkx),idslp(numsn)
      double precision histry(nhist,lastep+1)
c
c...  local variables
c
      integer i,idof,ilock,j,mode,ihist,n
c
cdebug      write(6,*) "Hello from printl_f!"
c
      if(numsn.eq.0) return
      do i=1,numsn
        n=idslp(i)
        idof=0
        ilock=1
        do j=1,ndof
          if(idx(j,n).ne.0) idof=idx(j,n)
        end do
        if(idof.eq.0) then
          ilock=1
        else if(nwinkx.ne.0) then
          do j=1,nwinkx
            if(iwinkx(2,j).eq.idof) then
              mode=iwinkx(1,j)
              if(mode.eq.0) ilock=0
              if(mode.eq.1) ilock=1
              if(mode.lt.0) then
                ihist=-mode
                ilock=0
                if(histry(ihist,nstep+1).ne.zero) ilock=1
              end if
            end if
          end do
        else
          ilock=0
        end if
        if(idsk.eq.2) write(kp) ilock
      end do
      return
      end
c
c version
c $Id: printl.f,v 1.3 2004/08/25 01:12:48 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
