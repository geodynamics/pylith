c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams
c  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
