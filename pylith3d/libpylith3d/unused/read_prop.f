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
      subroutine read_prop(prop,grav,dunits,yunits,viscunits,cohunits,
     & ndof,nprop,numat,idout,idsk,kr,kw,kp,ierr,ivisc,iplas,
     & prfile,ofile,pfile)
c
c...  subroutine to read data on element properties and
c     print a summary of the element group data.
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
      integer ndof,nprop,numat,idout,idsk,kr,kw,kp,ierr,ivisc,iplas
      double precision prop(nprop,numat),grav(3)
      double precision dunits,yunits,viscunits,cohunits
      character prfile*(*),ofile*(*),pfile*(*)
c
c...  included dimension and type statements
c
      include "labelp_dim.inc"
      include "labelc_dim.inc"
c
c...  intrinsic functions
c
      intrinsic index
c
c...  local variables
c
      integer i,j,inumat,n,nblks,nlast,ib,iln,iend
      character dummy*80
c
c...  included variable definitions
c
      include "labelp_def.inc"
      include "labelc_def.inc"
c
c...  open output files
c
      if(idout.gt.0) open(kw,file=ofile,status="old",access="append")
      if(idsk.eq.0) open(kp,file=pfile,status="old",access="append")
      if(idsk.eq.1) open(kp,file=pfile,status="old",access="append",
     & form="unformatted")
c
c...  open input file
c
      ierr=0
      if(idout.gt.0) write(kw,2010) numat
      open(kr,file=prfile,status="old",err=20)
c
c...  skip over units definitions.  Note that units must be specified
c     for density, young's modulus, viscosity coefficient, and cohesion,
c     even if they are not pertinent to the model.
c
      do i=1,4
        call pskip(kr)
        read(kr,"(a80)") dummy
        j=index(dummy,"=")
        if(j.eq.0) then
          ierr=2
          return
        end if
      end do
c
c...  read material properties
c
      call pskip(kr)
      do inumat=1,numat
        read(kr,*,err=30,end=30) n,(prop(i,n),i=1,nprop)
        prop(1,n)=yunits*prop(1,n)
        prop(3,n)=dunits*prop(3,n)
        prop(4,n)=viscunits*prop(4,n)
        prop(7,n)=cohunits*prop(7,n)
      end do
      close(kr)
c
c...  output plot results, if desired
c
      if(idsk.eq.0) then
        write(kp,1000) ivisc,iplas
        write(kp,1000) numat,nprop
        do n=1,numat
          do i=1,nprop
            write(kp,1010) prop(i,n)
          end do
        end do
      else if(idsk.eq.1) then
        write(kp) ivisc,iplas
        write(kp) numat,nprop
        write(kp) prop
      end if
c
c...  output ascii results, if desired
c
      if(idout.gt.0) then
        nblks=(numat-1)/5+1
        nlast=numat-5*(nblks-1)
        if(ivisc.eq.0) write(kw,2000)
        if(ivisc.eq.1) write(kw,2001)
        if(iplas.eq.0) write(kw,2002)
        if(iplas.eq.1) write(kw,2003)
        do ib=1,nblks
          write(kw,2020)
          iend=5
          if(ib.eq.nblks) iend=nlast
          write(kw,2023) (i+5*(ib-1),i=1,iend)
          do iln=1,nprop
            write(kw,2025) labelp(iln),(prop(iln,i+5*(ib-1)),i=1,iend)
          end do
        end do
        write(kw,2024)
        write(kw,2040)
        do i=1,ndof
          write(kw,2050) labelc(i),grav(i)
        end do
      end if
      if(idout.gt.0) close(kw)
      close(kp)
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
 1000 format(2i5)
 1010 format(e15.6)
 2010 format(1x///,
     x' m a t e r i a l   s e t   d a t a                    '  //5x,
     x' number of material sets . . . . . . . . . . (numat) =',i5)
 2000 format(//,' Viscous solution will not be used',//)
 2001 format(//,' Viscous solution will be used',//)
 2002 format(//,' Plastic solution will not be used',//)
 2003 format(//,' Plastic solution will be used',//)
 2020 format(//1x,' property ',31x,'material #')
 2023 format(16x,i5,4(9x,i5))
 2024 format(//,7x,' negative properties follow the given time history')
 2025 format(1x,a10,5(4x,1pe10.3))
 2040 format(1x//,
     x' a c c e l e r a t i o n   o f   g r a v i t y         ',//)
 2050 format(5x,
     x a4,'-direction  . . . . . . . . . . . . . . . . =',1pe17.5/)
c
      end
c
c version
c $Id: read_prop.f,v 1.1 2004/07/12 18:17:42 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
