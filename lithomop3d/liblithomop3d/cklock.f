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
      subroutine cklock(iwink,histry,ltim,lngth,nstep,nhist,lastep,
     & unlck)
c
c...checks whether fault locking forces are being imposed or removed
c   in the current time step.  if so, rff is forced by setting
c   ltim=.true.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      double precision cutoff
      parameter(cutoff=0.8d0)
c
c...  subroutine arguments
c
      integer lngth,nstep,nhist,lastep
      integer iwink(2,lngth)
      double precision histry(nhist,lastep+1)
      logical ltim,unlck
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer k,mode,ihist
      double precision diff
      logical test
c
cdebug      write(6,*) "Hello from cklock_f!"
c
      unlck=.false.
      test=.false.
      do k=1,lngth
        mode=iwink(1,k)
        if(mode.lt.0) then
          ihist=-mode
          diff=histry(ihist,nstep+1)-histry(ihist,nstep)
          if(abs(diff).gt.cutoff) ltim=.true.
          if(abs(diff).gt.cutoff) test=.true.
          if(diff.lt.-cutoff) unlck=.true.
        end if
        if(test.and.unlck) go to 10
      end do
10    continue
      return
      end
c
c version
c $Id: cklock.f,v 1.2 2004/06/21 19:09:38 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
