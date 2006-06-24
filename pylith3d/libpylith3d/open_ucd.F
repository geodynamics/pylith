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
      subroutine open_ucd(kucd,iprestress,nstep,ucdroot,iopt,iucd)
c
c...  Routine to determine filename for UCD output and open the file.
c
#include "binio.h"
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer kucd,iprestress,nstep,iopt,iucd
      character ucdroot*(*)
c
c...  external functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer i1,i2
      character filenm*200,cstep*5
c
cdebug      write(6,*) "Hello from open_ucd_f!"
c
      i1=nnblnk(ucdroot)
      i2=nchar(ucdroot)
      if(iopt.eq.1) then
        filenm=ucdroot(i1:i2)//".mesh.inp"
      else if(iopt.eq.2) then
        filenm=ucdroot(i1:i2)//".gmesh.inp"
      else if(iopt.eq.3) then
        if(iprestress.eq.izero) then
          write(cstep,"(i5.5)") nstep
        else
          cstep="prest"
        end if
        filenm=ucdroot(i1:i2)//".mesh.time."//cstep//".inp"
      else if(iopt.eq.4) then
        if(iprestress.eq.izero) then
          write(cstep,"(i5.5)") nstep
        else
          cstep="prest"
        end if
        filenm=ucdroot(i1:i2)//".gmesh.time."//cstep//".inp"
      end if
      if(iucd.eq.ione) then
        open(kucd,file=filenm,status="replace")
      else if(iucd.eq.itwo) then
        open(kucd,file=filenm,status="replace",
     &   BINARY_OPEN_SPEC,
     &   form="unformatted")
      end if
      return
      end
c
c version
c $Id: open_ucd.F,v 1.1 2005/08/05 19:58:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
