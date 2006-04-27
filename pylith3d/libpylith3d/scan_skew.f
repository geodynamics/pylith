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
      subroutine scan_skew(numrot,kr,rotation_units,skfile,ierr,
     & errstrng)
c
c...  subroutine to perform an initial scan of the skew rotation file
c     to determine the number of predefined skew rotations and the
c     units being used for rotations.
c
c     Note that the total number of skew rotations could be greater than
c     the value determined here if automatic rotation computation is being
c     performed for slippery nodes.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (this should not produce an
c             exception since skew BC are optional)
c         3:  Read error
c         5:  Units not specified
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
      integer numrot,kr,ierr
      character rotation_units*(*),skfile*(*),errstrng*(*)
c
c...  local constants
c
      character def(1)*14
      data def/"rotation_units"/
c
c...  local variables
c
      integer nget,j,n
      double precision skew(2)
      character units(1)*80
      logical units_defined(1)
c
c...  open input file
c
      ierr=izero
      numrot=izero
      nget=ione
      open(kr,file=skfile,status="old",err=10)
c
c...  get units, returning error 5 if they aren't found.
c
      call get_units(kr,nget,units_defined,units,def,ierr,errstrng)
      if(ierr.ne.izero) return
      rotation_units=units(1)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) n,(skew(j),j=1,nskdim)
        numrot=numrot+ione
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
        errstrng="scan_skew"
        close(kr)
        return
c
      end
c
c version
c $Id: scan_skew.f,v 1.2 2004/07/12 20:08:13 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
