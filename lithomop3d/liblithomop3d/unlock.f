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
      subroutine unlock(bwink,bintern,iwink,histry,nstep,nwink,nhist,
     & neq,lastep)
c
c       program to remove the winkler forces at a specified step.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nstep,nwink,nhist,neq,lastep
      integer iwink(2,nwink)
      double precision bwink(neq),bintern(neq),histry(nhist,lastep+1)
c
c...  local constants
c
      double precision cutoff
      parameter(cutoff=0.8d0)
c
c...  local variables
c
      integer k,idk,mode,ihist
      double precision diff
c
cdebug      write(6,*) "Hello from unlock_f!"
c
      do k=1,nwink
        mode=iwink(1,k)
        idk=iwink(2,k)
        if(mode.lt.izero) then
          ihist=-mode
          diff=histry(ihist,nstep+1)-histry(ihist,nstep)
          if(diff.lt.-cutoff) then
            bintern(idk)=bintern(idk)-bwink(idk)
            bwink(idk)=zero
          end if
        end if
      end do
      return
      end
c
c version
c $Id: unlock.f,v 1.4 2005/05/03 18:44:14 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
