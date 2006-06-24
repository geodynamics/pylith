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
      function rtsafe(x1,x2,xacc,ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,
     & efstsi,gam,dlam,deltp,alfap,iopd,idout,kto,kw,plas)
c
c...function to find the root of the effective stress function
c   adapted from Numerical Recipes
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iopd,idout,kto,kw
      double precision rtsafe,x1,x2,xacc,ae,bs,c,ds,dl1,dl2,t1,t2,emhu
      double precision anpwr,efstsi,gam,dlam,deltp,alfap
      logical plas
c
c...  defined constants
c
      include "rconsts.inc"
c
      integer maxit
      double precision eps
      parameter (maxit=100,eps=1.0d-18)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer j
      double precision xfac,acc,fl,df,fh,xl,xh,dxold,dx,f,temp
c
c*      write(6,*) "Hello from rtsafe_f!"
c
      xfac=half*(x1+x2)
      if(xfac.eq.zero) xfac=one
      acc=xacc*xfac
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x1,efstsi,
     & gam,dlam,fl,df,deltp,alfap,iopd,plas)
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,x2,efstsi,
     & gam,dlam,fh,df,deltp,alfap,iopd,plas)
      if((fl.gt.zero.and.fh.gt.zero).or.(fl.lt.zero.and.fh.lt.zero))
     & pause 'Root must be bracketed in rtsafe!'
      if(fl.eq.zero) then
        rtsafe=x1
        return
      else if(fh.eq.zero) then
        rtsafe=x2
        return
      else if(fl.lt.zero)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      end if
      rtsafe=half*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,rtsafe,efstsi,
     & gam,dlam,f,df,deltp,alfap,iopd,plas)
      do 11 j=1,maxit
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.zero
     *      .or. abs(two*f).gt.abs(dxold*df) ) then
          dxold=dx
          dx=half*(xh-xl)
          rtsafe=xl+dx
          if(abs(xl-rtsafe).le.eps)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(abs(temp-rtsafe).le.eps)return
        end if
        if(abs(dx).lt.acc) return
        call esf(ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,rtsafe,efstsi,
     &   gam,dlam,f,df,deltp,alfap,iopd,plas)
        if(f.lt.zero) then
          xl=rtsafe
        else if(f.gt.zero) then
          xh=rtsafe
        else
          return
        end if
11    continue
      write(kto,*) 'rtsafe exceeding maximum iterations!'
      if(idout.gt.1) write(kw,*) 'rtsafe exceeding maximum iterations!'
      return
      end
c
c version
c $Id: rtsafe.f,v 1.1 2004/07/12 19:35:27 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
