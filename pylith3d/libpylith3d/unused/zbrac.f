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
      subroutine zbrac(x1,x2,succes,ae,bs,c,ds,dl1,dl2,t1,t2,emhu,
     & anpwr,efstsi,gam,dlam,deltp,alfap,iopd,plas)
c
c...subroutine to bracket the root of the effective stress function
c     adapted from Numerical Recipes
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iopd
      double precision x1,x2,ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,efstsi
      double precision gam,dlam,deltp,alfap
      logical succes,plas
c
c...  defined constants
c
      include "rconsts.inc"
c
      integer ntry
      double precision factor
      parameter (factor=1.6d0,ntry=100)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer j
      double precision f1,df1,f2,df2
c
cdebug2      write(6,*) "Hello from zbrac_f!"
c
      if(x1.eq.x2)pause 'you have to guess an initial range'
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x1,efstsi,gam,dlam,
     & f1,df1,deltp,alfap,iopd,plas)
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x2,efstsi,gam,dlam,
     & f2,df2,deltp,alfap,iopd,plas)
cdebug2      write(6,*) "Point 1 in zbrac"
      succes=.true.
      do 11 j=1,ntry
        if(f1*f2.lt.zero)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+factor*(x1-x2)
          if(x1.lt.zero) x1=zero
          call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x1,efstsi,gam,
     &     dlam,f1,df1,deltp,alfap,iopd,plas)
        else
          x2=x2+factor*(x2-x1)
          call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x2,efstsi,gam,
     &     dlam,f2,df2,deltp,alfap,iopd,plas)
        end if
11    continue
      succes=.false.
      return
      end
c
c version
c $Id: zbrac.f,v 1.1 2004/07/12 21:25:39 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
