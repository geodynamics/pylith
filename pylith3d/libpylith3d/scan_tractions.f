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
      subroutine scan_tractions(numtractions,nsnodesmax,kr,
     & traction_units,tfile,ierr,errstrng)
c
c...  subroutine to perform an initial scan of the traction BC input
c     file to determine the total number of entries.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file (no exception should be raised
c             in this case since a traction BC file is optional)
c         3:  Read error
c       106:  Illegal element type
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
      integer numtractions,nsnodesmax,kr,ierr
      character traction_units*(*),tfile*(*),errstrng*(*)
c
c...  local constants
c
      character def(1)*14
      data def/"traction_units"/
c
c...  local variables
c
      integer nget,i,numverts
      integer tractionverts(4)
      double precision tractionvals(3)
      character units(1)*80
      logical units_defined(1)
cdebug      integer idb
c
cdebug      write(6,*) "Hello from scan_tractions_f!"
c
c...  open input file
c
      ierr=izero
      numtractions=izero
      numverts=izero
      nget=ione
      open(kr,file=tfile,status="old",err=10)
c
c...  get traction units, returning error 5 if they aren't found
c
      call get_units(kr,nget,units_defined,units,def,ierr,errstrng)
      if(ierr.ne.izero) return
      traction_units=units(1)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
c
c...  for now, determine number of vertices per face by evaluating read
c     errors
c
      read(kr,*,end=10,err=31) (tractionverts(i),i=1,4),
     & (tractionvals(i),i=1,ndof)
      numverts=4
      go to 33
 31   continue
        backspace(kr)
        read(kr,*,end=10,err=32) (tractionverts(i),i=1,3),
     &   (tractionvals(i),i=1,ndof)
        numverts=3
        go to 33
 32   continue
        ierr=106
        errstrng="scan_tractions"
        return
 33   continue
        backspace(kr)
 40   continue
        read(kr,*,end=10,err=30) (tractionverts(i),i=1,numverts),
     &   (tractionvals(i),i=1,ndof)
        numtractions=numtractions+1
        go to 40
c
c...  normal return
c
 10   continue
        close(kr)
        return
c
c...  read error
c
 30   continue
        ierr=3
        close(kr)
        errstrng="scan_tractions"
        return
c
      end
c
c version
c $Id: scan_tractions.f,v 1.6 2005/04/01 23:24:41 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
