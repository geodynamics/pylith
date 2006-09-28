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
      subroutine convert_case(string,iopt)
c
c...  routine to convert a character string from upper to lower case,
c     or vice-versa, depending on iopt.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iopt
      character *(*) string
c
c...  defined constants
c
      character upper*26,lower*26
      data lower/"abcdefghijklmnopqrstuvwxyz"/
      data upper/"ABCDEFGHIJKLMNOPQRSTUVWXYZ"/
c
c...  intrinsic functions
c
      intrinsic len
c
c...  local variables
c
      integer i,j,lenstr
c
      lenstr=len(string)
      if(iopt.eq.1) then
        do i=1,lenstr
          do j=1,26
            if(string(i:i).eq.upper(j:j)) then
              string(i:i)=lower(j:j)
              go to 10
            end if
          end do
 10       continue
        end do
      else
        do i=1,lenstr
          do j=1,26
            if(string(i:i).eq.lower(j:j)) then
              string(i:i)=upper(j:j)
              go to 20
            end if
          end do
 20       continue
        end do
      end if
c
      return
      end
c
c version
c $Id: convert_case.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
