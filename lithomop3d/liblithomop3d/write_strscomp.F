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
      subroutine write_strscomp(stol,dtol,epert,kw,idout,ofile)
c
c...subroutine to write parameters controlling the stress integration
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      double precision stol,dtol,epert
      integer kw,idout
      character ofile*(*)
c
c...  echo input to output file
c
      if(idout.gt.0) then
c
c...  at present, it is assumed that all f77 compilers support
c     access="append" and all f95 compilers support position="append".
c     If exceptions are found we will have to use FORTRAN_COMPILER_NAME
c     to differentiate.
c
#if defined FORTRAN_TYPE_F77
        open(kw,file=ofile,status="old",access="append")
#else
        open(kw,file=ofile,status="old",position="append")
#endif
        write(kw,700) stol,dtol,epert
        close(kw)
      end if
c
700   format(///' s t r e s s   i n t e g r a t i o n   ',
     & 'i n f o r m a t i o n'//,5x,
     & ' stress integration tolerance . . . . . . . (stol)  =',1pe15.8,
     & /,5x,
     & ' stress derivative tolerance. . . . . . . . (dtol)  =',1pe15.8,
     & /,5x,
     & ' maximum strain perturbation when computing',/,5x,
     & ' stress derivative. . . . . . . . . . . . . (epert) =',1pe15.8)
      return
      end
c
c version
c $Id: write_strscomp.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
