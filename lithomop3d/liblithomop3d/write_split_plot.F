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
      subroutine write_split_plot(idftn,numflt,kp,idsk,pfile)
c
c...  subroutine to write out each split node.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer numflt,kp,idsk
      integer idftn(numflt)
      character pfile*(*)
c
c...  local variables
c
      integer i
c
c...  write results to plot file, if requested
c
      if(idsk.eq.1) then
c
c...  at present, it is assumed that all f77 compilers support
c     access="append" and all f95 compilers support position="append".
c     If exceptions are found we will have to use FORTRAN_COMPILER_NAME
c     to differentiate.
c
#if defined FORTRAN_TYPE_F77
        open(kp,file=pfile,status="old",access="append")
#else
        open(kp,file=pfile,status="old",position="append")
#endif
        write(kp,"(i7)") numflt
        close(kp)
      else if(idsk.eq.2) then
#if defined FORTRAN_TYPE_F77
        open(kp,file=pfile,status="old",access="append",
     &   form="unformatted")
#else
        open(kp,file=pfile,status="old",position="append",
     &   form="unformatted")
#endif
        write(kp) numflt
        if(numflt.ne.0) write(kp) (idftn(i),i=1,numflt)
        close(kp)
      end if
c
      return
      end
c
c version
c $Id: write_split_plot.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
