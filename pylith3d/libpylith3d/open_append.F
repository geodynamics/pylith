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
      subroutine open_append(u,fin,sta,fm,ierr,errstrng,caller)
c
c...  subroutine to open an output file in append mode
c
c     The Autoconf 'configure' script determines whether the Fortran
c     compiler supports access="append", position="append", both, or
c     neither, at compile time.  Some compilers (e.g., gfortran) accept
c     both at compile time but reject one or the other at runtime;
c     therefore, this subroutine tries each available alternative in
c     turn.  If neither works, an error is returned.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer u,ierr
      character fin*(*),sta*(*),fm*(*),errstrng*(*),caller*(*)
c
c...  intrinsic functions
c
      intrinsic len
c
c...  local variables
c
      integer lenfm,jump
      save jump
      data jump/0/
c
      ierr=0
      lenfm=len(fm)
c
      goto (10,20) jump
c
 10   continue
#if defined FORTRAN_ACCESS_APPEND
      if(lenfm.eq.0) then
        open(u,file=fin,status=sta,access="append",err=20)
      else
        open(u,file=fin,status=sta,access="append",form=fm,err=20)
      end if
      jump=1
      goto 40
#endif
 20   continue
#if defined FORTRAN_POSITION_APPEND
      if(lenfm.eq.0) then
        open(u,file=fin,status=sta,position="append",err=30)
      else
        open(u,file=fin,status=sta,position="append",form=fm,err=30)
      end if
      jump=2
      goto 40
#endif
c
c...  error opening output file
c
 30   continue
        ierr=2
        errstrng=caller
        close(u)
        return
c
c...  normal return
c
 40   return
      end
c
c
c
      subroutine open_append_stop(u,fin,sta,fm,caller)
c
c...  convenience subroutine: open_append(), stop if error
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer u
      character fin*(*),sta*(*),fm*(*),caller*(*)
c
c...  local variables
c
      integer ierr
      character errstrng*256
c
      call open_append(u,fin,sta,fm,ierr,errstrng,caller)
      if(ierr.ne.0) then
        write(*,*) 'in function ', caller
        write(*,*) 'error opening file ', fin, ' for writing'
        stop
      end if
c
      return
      end
c
c End of file 
