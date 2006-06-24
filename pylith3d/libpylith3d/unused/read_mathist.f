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
      subroutine read_mathist(mhist,ivfamily,nvfamilies,npropsz,nhist,
     & kr,kw,kp,idout,idsk,mhfile,ofile,pfile,ierr,errstrng)
c
c...     reads material history information.
c
c     Error codes:
c         0:  No error
c         2:  Error opening output file
c         3:  Read error
c         4:  Write error
c       100:  Attempt to use undefined history
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nvfamilies,npropsz,nhist,kr,kw,kp,idout,idsk,ierr
      integer mhist(npropsz),ivfamily(5,nvfamilies)
      character mhfile*(*),ofile*(*),pfile*(*),errstrng*(*)
c
c...  local variables
c
      integer ifam,iprop,imhist,indprop,indpropg,i
c
      ierr=izero
      call ifill(mhist,izero,npropsz)
c
c...  open files
c
      open(kr,file=mhfile,status="old",err=10)
      if(idout.gt.izero) open(kw,file=ofile,err=40,status="old",
     & access="append")
      if(idsk.eq.ione) open(kp,file=pfile,err=40,status="old",
     & access="append")
      if(idsk.eq.itwo) open(kp,file=pfile,err=40,status="old",
     & access="append",form="unformatted")
      if(idout.gt.izero) write(kw,800)
c
c...  read material histories, if present
c
 70   continue
        call pskip(kr)
        read(kr,*,end=60,err=30) ifam,iprop,imhist
        if(imhist.le.izero.or.imhist.gt.nhist) then
          ierr=100
          errstrng="read_mathist"
          return
        end if
        indprop=ivfamily(5,ifam)
        indpropg=indprop+iprop-ione
        mhist(indpropg)=imhist
        if(idout.gt.izero) write(kw,810,err=50) ifam,iprop,imhist
        go to 70
 60   continue
      close(kr)
      if(idout.gt.izero) close(kw)
      if(idsk.eq.ione) write(kp,820,err=50) (mhist(i),i=1,npropsz)
      if(idsk.eq.itwo) write(kp,err=50) mhist
      close(kp)
c
c...  normal return
c
 10   continue
      return
c
c...  error reading input file
c
 30   continue
        ierr=3
        errstrng="read_mathist"
        close(kr)
        return
c
c...  error opening output file
c
 40   continue
        ierr=2
        errstrng="read_mathist"
        if(idout.gt.izero) close(kw)
        close(kp)
        return
c
c...  error writing to output file
c
 50   continue
        ierr=4
        errstrng="read_mathist"
        if(idout.gt.izero) close(kw)
        close(kp)
        return
c
800   format(//,
     & ' material property histories',/,
     & ' The specified properties for the following materials follow',/,
     & ' the given time history:',//,5x,
     & '      Material #     Property #     History #',/)
 810  format(6x,i5,8x,i5,8x,i5)
820   format(16i5)
c
      end
c
c version
c $Id: read_mathist.f,v 1.1 2005/04/20 19:00:14 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
