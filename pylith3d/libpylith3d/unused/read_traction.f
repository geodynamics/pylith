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
      subroutine read_traction(pres,pdir,trunits,ielno,iside,ihstry,
     & numpr,nen,npdir,kr,kw,idout,ierr,trfile,ofile)
c
c.......subroutine to read and print pressure load data
c
c        pres  = pressure array -- there is one entry per node on
c                the specified face
c        pdir  = pressure direction -- this has different meanings
c                depending on the geometry used:
c                plane strain/plane stress:  not used
c                axisymmetric:               not used
c                out-of-plane:               determines whether traction
c                                            is in-plane or out-of-plane
c                3D:                         specifies direction cosines
c                                            of traction
c        ielno = element number
c        iside = element face/side number
c        iphist = pressure load history #
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer numpr,nen,npdir,kr,kw,idout,ierr
      integer ielno(numpr),iside(numpr),ihstry(numpr)
      double precision pres(nen/2,numpr),pdir(npdir,numpr)
      double precision trunits
      character trfile*(*),ofile*(*)
c
c...  intrinsic functions
c
      intrinsic mod,index
c
c...  local variables
c
      integer i,n,npage
      character dummy*80
c
c...  open input file and skip past units specification.
c
      ierr=0
      if(numpr.eq.0) return
      open(kr,file=trfile,status="old",err=20)
      call pskip(kr)
      read(kr,"(a80)") dummy
      i=index(dummy,"=")
      if(i.eq.0) then
        ierr=2
        return
      end if
      call pskip(kr)
c
c...  read traction BC
c
      do n=1,numpr
        read(kr,*,end=30,err=30) ielno(n),iside(n),ihstry(n),
     &   (pres(i,n),i=1,nen/2),(pdir(i,n),i=1,npdir)
        do i=1,nen/2
          pres(i,n)=trunits*pres(i,n)
        end do
      end do
      close(kr)
c
c...  output results to ascii file, if requested
c
      if(idout.gt.0) then
        open(kw,file=ofile,status="old",access="append")
        npage=50
        do n=1,numpr
          if(n.eq.1.or.mod(n,npage).eq.0) write(kw,1000)
          write(kw,2000) ielno(n),iside(n),ihstry(n),
     &     (pres(i,n),i=1,nen/2),(pdir(i,n),i=1,npdir)
        end do
        close(kw)
      end if
c
c...  normal return
c
      return
c
c...  error opening input file
c
 20   continue
        ierr=1
        close(kr)
        return
c
c...  error reading input file
c
 30   continue
        ierr=3
        close(kr)
        return
c
1000  format(1x///,
     x' p r e s s u r e   l o a d   d a t a',//5x,
     x' element    side      hfac',8x,'pressure 1',5x,'pressure 2',/5x,
     x'  number   number',/)
 2000 format(6x,i7,4x,i7,4x,i7,7(2x,1pe15.8))
      end
c
c version
c $Id: read_traction.f,v 1.1 2004/07/12 19:14:38 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
