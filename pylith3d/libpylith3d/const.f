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
      subroutine const(maxstp,delt,alfa,maxit,ntdinit,lgdef,utol,
     & ftol,etol,itmax,nintg,i,naxstp,nfirst,rtimdat,deltp,alfap,
     & ntimdat,nstep,maxitp,ntdinitp,lgdefp,itmaxp,gtol)
c
c...this subroutine defines the constants for each time step group.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nintg,i,naxstp,nfirst
      integer maxstp(nintg),maxit(nintg),ntdinit(nintg)
      integer lgdef(nintg),itmax(nintg)
      double precision delt(nintg),alfa(nintg),utol(nintg),ftol(nintg)
      double precision etol(nintg)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "ntimdat_dim.inc"
      include "gtol_dim.inc"
c
c...  local variables
c
      integer ii
c
cdebug      write(6,*) "Hello from const_f!"
c
      ii=i
      deltp=delt(ii)
      rtimdat(1)=deltp
      alfap=alfa(ii)
      rtimdat(2)=alfap
      ntimdat(1)=nstep
      maxitp=maxit(ii)
      ntimdat(2)=maxitp
      ntdinitp=ntdinit(ii)
      ntimdat(3)=ntdinitp
      lgdefp=lgdef(ii)
      ntimdat(4)=lgdefp
      itmaxp=itmax(ii)
      ntimdat(5)=itmaxp
      gtol(1)=utol(ii)
      gtol(2)=ftol(ii)
      gtol(3)=etol(ii)
      nfirst=naxstp+1
      naxstp=nfirst+maxstp(ii)-1
      return
      end
c
c version
c $Id: const.f,v 1.3 2004/07/21 15:29:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
