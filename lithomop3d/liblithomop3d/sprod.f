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
      function sprod(r1,r2)
c
c...function to compute the scalar product of two 3x3 tensors (divided
c   by two).  The tensors are assumed to be symmetric, and thus are
c   given in vector form in the following order:  r(1,1), r(2,2),
c   r(3,3), r(1,2), r(2,3), and r(1,3).
c
      include "implicit.inc"
c
c...  function arguments
c
      double precision sprod,r1(6),r2(6)
c
c...  defined constants
c
      include "rconsts.inc"
c
c*      write(6,*) "Hello from sprod_f!"
c
      sprod=half*(r1(1)*r2(1)+r1(2)*r2(2)+r1(3)*r2(3))+r1(4)*r2(4)+
     & r1(5)*r2(5)+r1(6)*r2(6)
      return
      end
c
c version
c $Id: sprod.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
