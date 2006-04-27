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
      subroutine write_slip(nslip,numslp,numsn,kw,kp,idout,idsk,
     & ofile,pfile,ierr,errstrng)
c
c...  prints data on free slip interfaces
c
c     Error codes:
c         0:  No error
c         2:  Error opening output file
c         4:  Write error
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer numslp,numsn,kw,kp,idout,idsk,ierr
      integer nslip(nsdim,numslp)
      character ofile*(*),pfile*(*),errstrng*(*)
c
c...  included dimension and type statements
c
      include "labeld_dim.inc"
c
c...  intrinsic functions
      intrinsic mod
c
c...  local variables
c
      integer npage,i,j
c
c...  included variable definitions
c
      include "labeld_def.inc"
c
c...  open input file
c
      ierr=0
c
c...  open output files
c
      if(idout.gt.izero) then
        call open_append(kw,ofile,"old","",ierr,errstrng,"write_slip")
        if(ierr.ne.izero) return
      end if
      if(idsk.eq.ione) then
        call open_append(kp,pfile,"old","",ierr,errstrng,"write_slip")
        if(ierr.ne.izero) return
      else if(idsk.eq.itwo) then
        call open_append(kp,pfile,"old","unformatted",
     &                   ierr,errstrng,"write_slip")
        if(ierr.ne.izero) return
      end if
c
c...  loop over slippery node entries, and output results if desired
c
      if(idsk.eq.ione) write(kp,3000,err=60) numslp
      if(idsk.eq.itwo) write(kp,err=60) numslp
      if(numslp.ne.izero) then
        npage=50
        do i=1,numslp
          if((i.eq.ione.or.mod(i,npage).eq.izero).and.
     &     idout.gt.izero) then
            write(kw,1000,err=60) (labeld(j),j=1,ndof)
            write(kw,*,err=60) ' '
          end if
          if(idout.gt.izero) write(kw,2000,err=60) (nslip(j,i),
     &     j=1,ndof+2)
          if(idsk.eq.ione) write(kp,3000,err=60) (nslip(j,i),
     &     j=1,ndof+2)
        end do
      end if
      if(idsk.eq.itwo.and.numslp.ne.izero) write(kp) nslip
      if(idsk.eq.ione) write(kp,"(i7)",err=60) numsn
      if(idsk.eq.itwo) write(kp,err=60) numsn
      if(idout.gt.izero) close(kw)
      if(idsk.gt.izero) close(kp)
c
c...  normal return
c
      return
c
c...  error writing to output file
c
 60   continue
        ierr=4
        errstrng="write_slip"
        if(idout.gt.izero) close(kw)
        if(idsk.gt.izero) close(kp)
        return
c
 1000 format(//
     & 'i n t e r n a l   f r e e    s l i p  i n t e r f a c e   d a t
     &a'//
     & '  elem',4x,' node',12x,'dof weight'/
     & '   num',4x,'  num ',6(5x,a4))
 2000 format(i7,4x,i7,1x,6(5x,i4))
 3000 format(16i7)
      end
c
c version
c $Id: write_slip.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
