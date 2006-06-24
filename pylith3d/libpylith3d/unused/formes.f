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
      subroutine formes(x,d,dx,tfault,dmat,stn,skew,s,stemp,prop,gauss,
     & ien,lmx,lmf,iddmat,infin,n,ngauss,nddmat,nprop,nsd,ndof,nstr,nen,
     & nee,numnp,numfn,numslp,numrot,nskdim,iopt,ibbar,lgdef,idout,kto,
     & kw)
c
c...  subroutine to form the elemental stiffness matrix
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer infin,n,ngauss,nddmat,nprop,nsd,ndof,nstr,nen,nee,numnp
      integer numfn,numslp,numrot,nskdim,iopt,ibbar,lgdef,idout,kto,kw
      integer ien(nen),lmx(ndof,nen),lmf(nen)
      double precision x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn),dmat(nddmat,ngauss)
      double precision stn(nstr,ngauss),skew(nskdim,numnp),s(nee*nee)
      double precision stemp(nee*nee),prop(nprop),gauss(nsd+1,ngauss)
c
c...  included dimension and typing statements
c
      include "iddmat_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  local variables
c
cdebug      integer idb,jdb
      double precision xl(24),dl(24)
c
cdebug      write(6,*) "Hello from formes_f!"
cdebug      write(6,*) "From formes_f, ngauss,nsd,gauss: ",ngauss,nsd,
cdebug     & ((gauss(idb,jdb),idb=1,nsd+1),jdb=1,ngauss)
c
      call fill(s,zero,nee*nee)
c
c...  localize coordinates and update them for large deformations, if
c     required
c
      call lcoord(x,xl,ien,nen,nsd,numnp)
      if(lgdef.gt.izero) call ldupdat(d,dx,tfault,dl,xl,ien,lmx,lmf,
     & ndof,nsd,nen,numnp,numfn,numslp,iopt,lgdef)
c
c...  construct local stiffness matrix, symmetrize it, and rotate for
c     skew boundary conditions
c
      call stiffql(s,stn,dmat,xl,prop,gauss,ien,iddmat,infin,n,ngauss,
     & nddmat,nprop,ndof,nsd,nstr,nen,nee,ibbar,lgdef,idout,kto,kw)
      if(numrot.ne.izero) call rstiff(s,stemp,skew,ien,ndof,numnp,nen,
     & nee,nskdim)
      return
      end
c
c version
c $Id: formes.f,v 1.1 2004/06/15 20:03:35 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
