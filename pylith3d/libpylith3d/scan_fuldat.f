c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
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
      subroutine scan_fuldat(icode,lastep,icontr,kr,fofile,ierr,
     & errstrng)
c
c...  subroutine to perform an initial scan of the file specifying the
c     timesteps for which to produce output to determine the number of
c     outputs to produce.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (an exception should be raised if
c             a full solution is being performed and the number of time
c             steps is greater than zero)
c         3:  Read error
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer icode,lastep,icontr,kr,ierr
      character fofile*(*),errstrng*(*)
c
c...  local variables
c
      integer iprint
c
c...  open input file
c
      ierr=izero
      icontr=izero
      if(lastep.eq.izero) return
      open(kr,file=fofile,status="old",err=20)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) iprint
        icontr=icontr+ione
        go to 40
c
c...  normal return
c
 10   continue
        if(icontr.eq.izero.and.icode.eq.ithree) then
          ierr=3
          errstrng="scan_fuldat"
        end if
        close(kr)
        return
c
c...  error opening file
c
 20   continue
	if(icode.eq.ithree) then
          ierr=1
          errstrng="scan_fuldat"
        end if
        close(kr)
        return
c
c...  read error
c
 30   continue
        ierr=3
        errstrng="scan_fuldat"
        close(kr)
        return
c
      end
c
c version
c $Id: scan_fuldat.f,v 1.2 2004/07/12 19:58:02 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
