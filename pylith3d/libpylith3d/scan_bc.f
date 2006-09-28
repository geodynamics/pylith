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
      subroutine scan_bc(numbc,kr,displacement_units,
     & velocity_units,force_units,bcfile,ierr,errstrng)
c
c...  subroutine to perform an initial scan of the boundary condition
c     file to determine the number of boundary condition entries and the
c     units being used for displacement, velocity, and force.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
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
      integer numbc,kr,ierr
      character displacement_units*(*),velocity_units*(*)
      character force_units*(*),bcfile*(*),errstrng*(*),errtmp*512
c
c...  local constants
c
      character def(3)*18
      data def/"displacement_units","velocity_units","force_units"/
c
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer nget,j,n,i1,i2
      integer ibond(3)
      double precision bond(3)
      character units(3)*80
      logical units_defined(3)
c
c...  open input file
c
      ierr=izero
      numbc=izero
      nget=ithree
      open(kr,file=bcfile,status="old",err=20)
c
c...  get units, returning error 2 if they aren't found.
c
      call get_units(kr,nget,units_defined,units,def,ierr,errtmp)
      if(ierr.ne.izero) then
        i1=nnblnk(errtmp)
        i2=nchar(errtmp)
        errstrng="scan_bc:"//errtmp(i1:i2)
        return
      end if
      displacement_units=units(1)
      velocity_units=units(2)
      force_units=units(3)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) n,(ibond(j),j=1,ndof),
     &   (bond(j),j=1,ndof)
        numbc=numbc+1
        go to 40
c
c...  normal return
c
 10   continue
        close(kr)
        return
c
c...  error opening file
c
 20   continue
        close(kr)
	ierr=1
        errstrng="scan_bc"
        return
c
c...  read error
c
 30   continue
        close(kr)
        ierr=3
        errstrng="scan_bc"
        return
c
      end
c
c version
c $Id: scan_bc.f,v 1.3 2004/08/12 20:49:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
