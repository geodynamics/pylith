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
      subroutine scan_prop(nprop,numat,kr,ierr,density_units,
     & young_units,viscosity_coeff_units,cohesion_units,
     & ivisc,iplas,imhist,pfile)
c
c...  subroutine to perform an initial scan of the material properties
c     definition file to determine the number of material types and the
c     units being used for density, volume change, young's modulus,
c     viscosity coefficient, and cohesion.
c     this routine also sets the viscous and plastic solution flags.
c
c     Exceptions should be produced for all error codes since properties
c     must be specified.
c
c     Error codes:
c         0:  No error
c         1:  Error opening input file
c         2:  Units not specified
c         3:  Read error
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nprop,numat,kr,ierr,ivisc,iplas,imhist
      character density_units*(*),young_units*(*)
      character viscosity_coeff_units*(*),cohesion_units*(*),pfile*(*)
c
c...  defined constants
c
      character def(4)*27
      data def/"density_units","youngs_modulus_units",
     & "viscosity_coefficient_units","cohesion_units"/
c
c...  local variables
c
      integer nget,j,n
      double precision prop(20)
      character units(4)*80
      logical units_defined(4)
c
c...  open input file
c
      ierr=0
      numat=0
      nget=4
      open(kr,file=pfile,status="old",err=20)
c
c...  get units, returning error 2 if they aren't found.
c
      call get_units(kr,ierr,nget,units_defined,units,def)
      if(ierr.eq.2) return
      density_units=units(1)
      young_units=units(2)
      viscosity_coeff_units=units(3)
      cohesion_units=units(4)
c
c... scan the file, counting the number of entries.
c    Note:  Due to speed considerations, we are not allowing the option
c    of putting comments within the list.  To do this, we
c    would need to scan each line twice to determine whether it was a
c    comment.
c
      ivisc=0
      iplas=0
      imhist=0
      call pskip(kr)
 40   continue
        read(kr,*,end=10,err=30) n,(prop(j),j=1,nprop)
        call matflg(prop,nprop,ivisc,iplas,imhist)
        numat=numat+1
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
        return
c
c...  read error
c
 30   continue
        ierr=3
        close(kr)
        return
c
      end
c
c version
c $Id: scan_prop.f,v 1.1 2004/07/12 20:04:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
