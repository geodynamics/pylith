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
      subroutine write_strscomp(stol,dtol,epert,kw,idout,ofile)
c
c...subroutine to write parameters controlling the stress integration
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      double precision stol,dtol,epert
      integer kw,idout
      character ofile*(*)
c
c...  echo input to output file
c
      if(idout.gt.0) then
        call open_append_stop(kw,ofile,"old","","write_strscomp")
        write(kw,700) stol,dtol,epert
        close(kw)
      end if
c
700   format(///' s t r e s s   i n t e g r a t i o n   ',
     & 'i n f o r m a t i o n'//,5x,
     & ' stress integration tolerance . . . . . . . (stol)  =',1pe15.8,
     & /,5x,
     & ' stress derivative tolerance. . . . . . . . (dtol)  =',1pe15.8,
     & /,5x,
     & ' maximum strain perturbation when computing',/,5x,
     & ' stress derivative. . . . . . . . . . . . . (epert) =',1pe15.8)
      return
      end
c
c version
c $Id: write_strscomp.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
