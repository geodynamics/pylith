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
      subroutine get_units(ifile,nget,units_defined,units,def,ierr,
     & errstrng)
c
c...  routine to get the variable units for variable(s) 'def'.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ifile,nget,ierr
      character units(nget)*(*),def(nget)*(*),errstrng*(*)
      logical units_defined(nget)
c
c...  intrinsic functions
c
      intrinsic index
c
c...  user-defined functions
c
      integer nchar,nnblnk
      external nchar,nnblnk
c
c...  local variables
c
      integer lendef,iopt,i,ii,j,jj,ib,ie
      character string*80,substr1*80,substr2*80
c
      ierr=0
      iopt=1
      do i=1,nget
        units_defined(i)=.false.
      end do
c
      do i=1,nget
        call pskip(ifile)
        read(ifile,"(a80)") string
        ii=index(string,"=")
        if(ii.eq.0) then
          ierr=5
          errstrng="get_units"
          return
        end if
        substr1=string(1:ii-1)
        call convert_case(substr1,iopt)
        do j=1,nget
          lendef=nchar(def(j))
          jj=index(substr1,def(j)(1:lendef))
          if(jj.ne.0) then
            units_defined(j)=.true.
            substr2=string(ii+1:)
            ib=nnblnk(substr2)
            ie=nchar(substr2)
            if(ie.eq.0) then
              ierr=5
              errstrng="get_units"
              return
            end if
            units(j)=substr2(ib:ie)
          end if
        end do
      end do
c
      do i=1,nget
        if(.not.units_defined(i)) then
          ierr=5
          errstrng="get_units"
          return
        end if
      end do
c
      return
      end
c
c version
c $Id: get_units.f,v 1.2 2004/07/06 20:33:36 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
