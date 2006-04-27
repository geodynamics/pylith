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
      subroutine matprtb(stn,eps,beta,betb,dmat,prop,rtimdat,
     & rgiter,iddmat,n,nstr,ndof,nprop,nddmat,ipstrs,nstep,lgdefp,idout,
     & kto,kw,ivisc,iplas)
c
c...constructs the tangent constitutive matrix by perturbing the
c   strain vector about its current state.  The corresponding stress
c   perturbations are computed using the effective stress function
c   algorithm.  The numerical derivative algorithm is adapted from
c   Numerical Recipes.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer n,nstr,ndof,nprop,nddmat,ipstrs,nstep,lgdefp,idout,kto,kw
      integer ivisc,iplas
      double precision stn(nstr),eps(nstr),beta(nstr),betb(nstr)
      double precision dmat(nddmat),prop(nprop)
c
c...  included dimension and type statements
c
      include "iddmat_dim.inc"
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
      integer ntab
      double precision con,con2,big,safe
      parameter(con=1.4d0,con2=con*con,big=1.d30,ntab=10,safe=2.0d0)
c
c...  intrinsic functions
c
      intrinsic abs,max
c
c...  user-defined functions
c
      double precision sprod
      external sprod
c
c...  local variables
c
      integer i,j,ii,jj
      double precision pinit,hh,fac
      double precision a(6,ntab,ntab),dtmp(6),dbeta(6),dbetb(6),stnb(6)
      double precision stnf(6),ee(6),errt(6),err(6)
      logical*4 converge
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
c
c*      write(6,*) "Hello from matprtb_f!"
c
      call fill(dmat,zero,nddmat)
      pinit=epert*sprod(eps,eps)
      if(pinit.lt.dtol) pinit=epert
      do i=1,6
        hh=pinit
        call dcopy(6,eps,ione,ee,ione)
        ee(4)=half*ee(4)
        ee(5)=half*ee(5)
        ee(6)=half*ee(6)
        ee(i)=ee(i)-hh
        call esfcomp(stn,stnb,ee,beta,dbeta,betb,dbetb,prop,
     &   rtimdat,stol,n,nstr,ndof,nprop,ipstrs,nstep,lgdefp,idout,kto,
     &   kw,ivisc,iplas)
        ee(i)=ee(i)+two*hh
        call esfcomp(stn,stnf,ee,beta,dbeta,betb,dbetb,prop,
     &   rtimdat,stol,n,nstr,ndof,nprop,ipstrs,nstep,lgdefp,idout,kto,
     &   kw,ivisc,iplas)
        do j=1,i
          a(j,1,1)=(stnf(j)-stnb(j))/(two*hh)
          err(j)=big
        end do
        do ii=2,ntab
          hh=hh/con
          call dcopy(6,eps,ione,ee,ione)
          ee(4)=half*ee(4)
          ee(5)=half*ee(5)
          ee(6)=half*ee(6)
          ee(i)=ee(i)-hh
          call esfcomp(stn,stnb,ee,beta,dbeta,betb,dbetb,prop,
     &     rtimdat,stol,n,nstr,ndof,nprop,ipstrs,nstep,lgdefp,idout,kto,
     &     kw,ivisc,iplas)
          ee(i)=ee(i)+two*hh
          call esfcomp(stn,stnf,ee,beta,dbeta,betb,dbetb,prop,
     &     rtimdat,stol,n,nstr,ndof,nprop,ipstrs,nstep,lgdefp,idout,kto,
     &     kw,ivisc,iplas)
          converge=.true.
          do j=1,i
            a(j,1,ii)=(stnf(j)-stnb(j))/(two*hh)
            fac=con2
            do jj=2,ii
              a(j,jj,ii)=(a(j,jj-1,ii)*fac-a(j,jj-1,ii-1))/(fac-one)
              fac=con2*fac
              errt(j)=max(abs(a(j,jj,ii)-a(j,jj-1,ii)),abs(a(j,jj,ii)-
     &         a(j,jj-1,ii-1)))
              if(errt(j).le.err(j)) then
                err(j)=errt(j)
                dtmp(j)=a(j,jj,ii)
              end if
            end do
            converge=converge.and.(abs(a(j,ii,ii)-a(j,ii-1,ii-1)).ge.
     &       safe*err(j))
          end do
          if(converge) go to 10
        end do
10      continue
        do j=1,i
          if(j.gt.3) dtmp(j)=half*dtmp(j)
          dmat(iddmat(j,i))=dtmp(j)
        end do
      end do
      return
      end
c
c version
c $Id: matprtb.f,v 1.1 2004/07/07 15:44:31 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
