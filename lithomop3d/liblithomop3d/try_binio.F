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
      subroutine try_binio(u,ierr,errstrng)
c
c...  subroutine to try binary I/O to see if it works
c
#include "binio.h"
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer u,ierr
      character errstrng*(*)
c
c...  local variables
c
      integer intlen,floatlen,p
      integer testint
      real testfloat
c
      ierr=0
c
      intlen=0
      floatlen=0
      p=0
      testint=42
      testfloat=4.2
c
#if defined FORTRAN_NO_BINARY
      ierr=400
      goto 10
#else
c
c...  open a scratch file
c
      open(u,status="scratch",
     & BINARY_OPEN_SPEC,
     & form="unformatted",err=402)
c
c...  write an integer and a float and measure their size
c
      write(u,BINARY_POS=1,BINARY_ERR(404)) testint
      inquire(u,BINARY_INQ_POS=p)
      intlen = p - 1
      write(u,BINARY_POS=p,BINARY_ERR(404)) testfloat
      inquire(u,BINARY_INQ_POS=p)
      floatlen = p - intlen - 1
      if(intlen.ne.ifour.or.floatlen.ne.ifour) then
        ierr=405
        goto 10
      end if
c
c...  read them back and verify that the values are the same
c
      testint=0
      testfloat=0.0
      read(u,BINARY_POS=1,BINARY_ERR(403)) testint
      p = intlen + 1
      read(u,BINARY_POS=p,BINARY_ERR(403)) testfloat
      close(u)
      if(testint.ne.42.or.testfloat.ne.4.2) then
        ierr=406
        goto 10
      end if
#endif
c
c...  normal return
c
      return
c
c...  error
c
 402  ierr=402
      goto 10
 403  ierr=403
      goto 10
 404  ierr=404
      goto 10
c
 10   errstrng="try_binio"
      close(u)
      return
c
      end
c
c End of file 
