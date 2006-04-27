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
      subroutine scan_split(numfn,kr,spfile,ierr,errstrng)
c
c...  subroutine to perform an initial scan of the split node input
c     file to determine the number of split node entries (numfn).  Note
c     that the number of actual split nodes (numflt) is determined in
c     the parsing routine.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (no exception should be raised
c             in this case since a split node file is optional)
c         3:  Read error
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
      integer numfn,kr,ierr
      character spfile*(*),errstrng*(*)
c
c...  local variables
c
      integer j
      integer nfault(3)
      double precision fault(3)
c
c...  open input file
c
      ierr=izero
      numfn=izero
      open(kr,file=spfile,status="old",err=10)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) (nfault(j),j=1,3),(fault(j),j=1,ndof)
        numfn=numfn+ione
        go to 40
c
c...  normal return
c
 10   continue
        close(kr)
        return
c
c...  error opening file -- not used since this file is optional
c
c* 20   continue
c*        ierr=1
c*        close(kr)
c*        return
c
c...  read error
c
 30   continue
        ierr=3
        errstrng="scan_split"
        close(kr)
        return
c
      end
c
c version
c $Id: scan_split.f,v 1.2 2004/07/12 20:12:38 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
