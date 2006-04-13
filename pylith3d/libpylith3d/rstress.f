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
      subroutine rstress(stn,r)
c
c...routine to rotate the stresses using the rotation matrix r
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      double precision stn(6),r(3,3)
c
c...  defined constants
c
      include "rconsts.inc"
      integer id
      parameter(id=3)
c
c...  local variables
c
      double precision stnm(3,3),tmp(3,3)
c
cdebug      write(6,*) "Hello from rstress_f!"
c
      stnm(1,1)=stn(1)
      stnm(2,2)=stn(2)
      stnm(3,3)=stn(3)
      stnm(1,2)=stn(4)
      stnm(2,3)=stn(5)
      stnm(1,3)=stn(6)
      stnm(2,1)=stnm(1,2)
      stnm(3,1)=stnm(1,3)
      stnm(3,2)=stnm(2,3)
      call dgemm("n","n",id,id,id,one,r,id,stnm,id,zero,tmp,id)
      call dgemm("n","t",id,id,id,one,tmp,id,r,id,zero,stnm,id)
      stn(1)=stnm(1,1)
      stn(2)=stnm(2,2)
      stn(3)=stnm(3,3)
      stn(4)=stnm(1,2)
      stn(5)=stnm(2,3)
      stn(6)=stnm(1,3)
      return
      end
c
c version
c $Id: rstress.f,v 1.2 2004/08/12 02:25:15 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
