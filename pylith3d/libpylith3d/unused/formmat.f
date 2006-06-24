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
      subroutine formmat(
     & stn,dmat,eps,beta,betb,                                          ! stress
     & prop,mat,iddmat,                                                 ! elemnt
     & histry,rtimdat,ntimdat,rgiter,iter,                              ! timdat
     & ngauss,nddmat,ndmat,nprop,numat,ndof,nstr,numel,ipstrs,nhist,    ! dimens
     & lastep,idout,kto,kw)                                             ! dimens
c
c...subroutine to reform the d-matrix
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iter,ngauss,nddmat,ndmat,nprop,numat,ndof,nstr,numel
      integer ipstrs,nhist,lastep,idout,kto,kw
      integer mat(numel)
      double precision stn(nstr,ngauss,numel),dmat(nddmat,ngauss,ndmat)
      double precision eps(nstr,ngauss,numel),beta(nstr,ngauss,numel)
      double precision betb(nstr,ngauss,numel),prop(nprop,numat)
      double precision histry(nhist,lastep+1)
c
c...  included dimension and type statements
c
      include "iddmat_dim.inc"
      include "ntimdat_dim.inc"
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
c
c...  local variables
c
cdebug      integer idb,jdb,kdb
      integer m,l,n
      double precision ptmp(30)
c
c...  included variable definitions
c
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "From formmat_f, ngauss: ",ngauss
cdebug      write(6,*) "Hello from formmat_f!"
cdebug2      write(6,*) "iter,ngauss,nddmat,ndmat,nprop,numat,ndof,nstr,numel"
cdebug2      write(6,*) "ipstrs,nhist,lastep,idout,kto,kw"
cdebug2      write(6,*) iter,ngauss,nddmat,ndmat,nprop,numat,ndof,nstr,numel
cdebug2      write(6,*) ipstrs,nhist,lastep,idout,kto,kw
cdebug2      write(6,*) "stn:",(((stn(idb,jdb,kdb),idb=1,nstr),jdb=1,ngauss),
cdebug2     & kdb=1,numel)
cdebug2      write(6,*) "dmat:",(((dmat(idb,jdb,kdb),idb=1,nddmat),
cdebug2     & jdb=1,ngauss),kdb=1,ndmat)
cdebug2      write(6,*) "eps:",(((eps(idb,jdb,kdb),idb=1,nstr),jdb=1,ngauss),
cdebug2     & kdb=1,numel)
cdebug2      write(6,*) "beta:",(((beta(idb,jdb,kdb),idb=1,nstr),jdb=1,ngauss),
cdebug2     & kdb=1,numel)
cdebug2      write(6,*) "prop:",((prop(idb,jdb),idb=1,nprop),jdb=1,numat)
cdebug2      write(6,*) "histry:",((histry(idb,jdb),idb=1,nhist),
cdebug2     & jdb=1,lastep+1)
c
      if((ivisc.eq.izero).and.(iplas.eq.izero).and.iter.eq.ione) then
        do m=1,numat
          call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &     idout,kto,kw,imhist)
          do l=1,ngauss
            call matinit(stn(1,l,m),eps(1,l,m),beta(1,l,m),betb(1,l,m),
     &       dmat(1,l,m),ptmp,rtimdat,iddmat,nstr,nddmat,nprop,
     &       ndof,ipstrs,nstep,lgdefp,ivisc,iplas)
cdebug            write(6,"(2i7,21(2x,1pe15.8))") m,l,(dmat(jdb,l,m),jdb=1,nddmat)
          end do
        end do
      else if(((ivisc.eq.ione).or.(iplas.eq.ione)).and.
     & (iter.eq.ione.or.(iter.gt.ione.and.nstep.gt.izero))) then
        do n=1,numel
          m=mat(n)
          call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &     idout,kto,kw,imhist)
          do l=1,ngauss
            if(iter.eq.ione) call matinit(stn(1,l,n),eps(1,l,n),
     &       beta(1,l,n),betb(1,l,n),dmat(1,l,n),ptmp,rtimdat,
     &       iddmat,nstr,nddmat,nprop,ndof,ipstrs,nstep,lgdefp,ivisc,
     &       iplas)
            if(iter.gt.ione) call matprtb(stn(1,l,n),eps(1,l,n),
     &       beta(1,l,n),betb(1,l,n),dmat(1,l,n),ptmp,rtimdat,
     &       rgiter,iddmat,n,nstr,ndof,nprop,nddmat,ipstrs,nstep,lgdefp,
     &       idout,kto,kw,ivisc,iplas)
cdebug2            write(6,*) "n,l,dmat:",n,l,(dmat(jdb,l,n),jdb=1,nddmat)
          end do
        end do
      end if
      return
      end
c
c version
c $Id: formmat.f,v 1.1 2004/07/02 18:42:02 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
