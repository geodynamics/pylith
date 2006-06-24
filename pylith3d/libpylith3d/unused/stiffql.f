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
      subroutine stiffql(s,stn,dmat,xl,prop,gauss,ien,iddmat,infin,n,
     & ngauss,nddmat,nprop,ndof,nsd,nstr,nen,nee,ibbarp,lgdefp,idout,
     & kto,kw)
c
c...computes the local stiffness matrix at 2 x 2 x 2 barlow points
c   k=(b)t*d*b.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer infin,n,ngauss,nddmat,nprop,ndof,nsd,nstr,nen,nee,ibbarp
      integer lgdefp,idout,kto,kw
      integer ien(nen)
      double precision s(nee,nee),stn(nstr,ngauss),dmat(nddmat,ngauss)
      double precision xl(nsd,nen),prop(nprop),gauss(nsd+1,ngauss)
c
c...  included dimension and type statements
c
      include "iddmat_dim.inc"
c
c...  defined constants
c
      include "rconsts.inc"
c
c...  local variables
c
cdebug      integer idb,jdb
      integer iopt,l,ibbart
      double precision ,vol
      double precision sh(4,8,8),shbar(4,8),b(6,24),db(6,24),det(8)
      double precision dtmp(6,6)
c
c...form shape functions for each integration point, then compute
c   average for b-bar formalism.
c
c*      write(6,*) "Hello from stiffql_f!"
cdebug      write(6,*) "From stiffql_f, ngauss,nsd,gauss: ",ngauss,nsd,
cdebug     & ((gauss(idb,jdb),idb=1,nsd+1),jdb=1,ngauss)
c
      vol=zero
      iopt=2
      do l=1,ngauss
        call shapql(gauss(1,l),xl,det(l),sh(1,1,l),ien,nen,nsd,infin,
     &   iopt,n,idout,kto,kw)
        det(l)=gauss(4,l)*det(l)
        vol=vol+det(l)
      end do
      ibbart=0
      if(ibbarp.eq.1.and.ngauss.gt.1) then
        call meanshql(shbar,sh,vol,det,ngauss)
        ibbart=1
      end if
c
c...construct b matrix, then form intermediate db=dmat*b, and finally
c   the stiffness (b)t*dmat*b multiplied by appropriate weight for
c   integral over element
c*  if lgdefp=2, compute additional contribution due to initial stress
c
      if(lgdefp.eq.2) call stiffldql(sh,stn,s,det,ngauss,ndof,nstr,nen,
     & nee)
      do l=1,ngauss
        call bmatrixql(b,sh(1,1,l),shbar,ibbart)
        call getmat(dtmp,dmat(1,l),iddmat,nstr,nddmat)
        call dsymm("l","l",nstr,nee,det(l),dtmp,nstr,b,nstr,zero,db,
     &   nstr)
        call dgemm("t","n",nee,nee,nstr,one,b,nstr,db,nstr,one,s,nee)
c*******  two things to do:  1.  See if I can use triangular dmat.
c*******                     2.  See if there is a LAPACK (or maybe BLAS
c*******                         routine) to perform symmetric
c*******                         matrix multiplication.
      end do
cdebug      write(6,*) "From stiffql_f, n,s: ",n,
cdebug     & ((s(idb,jdb),idb=1,nee),jdb=1,nee)
      return
      end
c
c version
c $Id: stiffql.f,v 1.1 2004/06/15 20:03:35 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
