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
      subroutine write_subiter(rmult,rmin,nsol,maxcg,kw,idout,ofile)
c
c...subroutine to write parameters controlling the subiteration
c   process
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nsol,maxcg,kw,idout
      include "rmult_dim.inc"
      include "rmin_dim.inc"
      character ofile*(*)
c
c...  local variables
c
      integer i
c
      if(idout.gt.0) then
        open(kw,file=ofile,status="old",access="append")
        write(kw,700) nsol,maxcg,(rmult(i),i=1,3),(rmin(i),i=1,3)
        close(kw)
      end if
c
      return
c
700   format(///' l i n e a r   s o l u t i o n   ',
     & 'i n f o r m a t i o n'//,5x,
     & ' PCG solution method  . . . . . . . . . . . (nsol)  =',i5,//,5x,
     & '    eq.1,  diagonal preconditioner,                  ',/,5x,
     & '           use zero initial guess                    ',/,5x,
     & '    eq.2,  diagonal preconditioner,                  ',/,5x,
     & '           use displacements from previous step      ',/,5x,
     & '    eq.3,  symmetrized Gauss-Seidel preconditioner   ',/,5x,
     & '           use zero initial guess                    ',/,5x,
     & '    eq.4,  symmetrized Gauss-Seidel preconditioner   ',/,5x,
     & '           use displacements from previous step      ',//,5x,
     & ' maximum number of PCG iterations . . . . . (maxcg) =',i5,//,5x,
     & ' accuracy multiplier for displacement . . . (rmult1)=',1pe15.8,
     & /,5x,
     & ' accuracy multiplier for force. . . . . . . (rmult2)=',1pe15.8,
     & /,5x,
     & ' accuracy multiplier for energy . . . . . . (rmult3)=',1pe15.8,
     & /,5x,
     & ' minimum displacement accuracy. . . . . . . (rmin1) =',1pe15.8,
     & /,5x,
     & ' minimum force accuracy . . . . . . . . . . (rmin2) =',1pe15.8,
     & /,5x,
     & ' minimum energy accuracy. . . . . . . . . . (rmin3) =',1pe15.8)
      end
