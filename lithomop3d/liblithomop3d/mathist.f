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
      subroutine mathist(ptmp,prop,histry,nprop,m,nstep,nhist,lastep,
     & idout,kto,kw,imhist)
c
c...subroutine to assign material properties based on time histories
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nprop,m,nstep,nhist,lastep,idout,kto,kw,imhist
      double precision ptmp(nprop),prop(nprop),histry(nhist,lastep+1)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  intrinsic functions
c
      intrinsic nint,abs
c
c...  local variables
c
      integer ihist,i
c
c*      write(6,*) "Hello from mathist_f!"
c
      if(imhist.eq.1) then
        ihist=nint(prop(nprop))
        if(ihist.gt.nhist.or.ihist.lt.0) then
          if(idout.gt.1) write(kw,1000) ihist,m
          write(kto,1000) ihist,m
          stop
        end if
        do i=1,nprop-1
          ptmp(i)=abs(prop(i))
          if(prop(i).lt.zero.and.ihist.ne.0)
     &     ptmp(i)=histry(ihist,nstep+1)*abs(prop(i))
        end do
      else
	call dcopy(nprop,prop,ione,ptmp,ione)
      end if
c
1000  format(//' fatal material property error!'//
     & ' attempt to use undefined load history # ',i5,
     & ' for material set # ',i5)
      return
      end
c
c version
c $Id: mathist.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
