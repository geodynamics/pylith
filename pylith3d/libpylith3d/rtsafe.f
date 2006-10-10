c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      FUNCTION rtsafe(func,x1,x2,xacc,rpar,nrpar,ipar,nipar,
     & ierr,errstrng)
c
c...  function to find the root of a real function.  Modified from
c     Numerical Recipes.
c
c     Error codes:
c         0:  No error
c       115:  Root not initially bracketed
c       411:  Maximum iterations exceeded
c
      include "implicit.inc"
c
c...  parameter definitions
c
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
c
c...  subroutine arguments
c
      integer nrpar,nipar,ierr
      integer ipar(nipar)
      double precision rtsafe,x1,x2,xacc,rpar(nrpar)
      character errstrng*(*)
c
c...  external routines
c
      EXTERNAL func
c
c...  local variables
c
      INTEGER j
      double precision df,dx,dxold,f,fh,fl,temp,xh,xl
c
      ierr=0
c
      call func(x1,fl,df,rpar,nrpar,ipar,nipar)
      call func(x2,fh,df,rpar,nrpar,ipar,nipar)
      if((fl.gt.0.0d0.and.fh.gt.0.0d0).or.
     & (fl.lt.0.0d0.and.fh.lt.0.0d0)) then
        ierr=115
        errstrng="rtsafe"
        return
      end if
c
      if(fl.eq.0.0d0)then
        rtsafe=x1
        return
      else if(fh.eq.0.0d0)then
        rtsafe=x2
        return
      else if(fl.lt.0.0d0)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
c
      rtsafe=0.5d0*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      call func(rtsafe,f,df,rpar,nrpar,ipar,nipar)
c
      do j=1,MAXIT
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0.0d0.or.
     &   abs(2.0d0*f).gt.abs(dxold*df) ) then
          dxold=dx
          dx=0.5d0*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.xacc) return
        call func(rtsafe,f,df,rpar,nrpar,ipar,nipar)
        if(f.lt.0.0d0) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
      end do
      ierr=411
      errstrng="rtsafe"
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .zW.
