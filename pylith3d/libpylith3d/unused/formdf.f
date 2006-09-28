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
      subroutine formdf(b,x,d,dx,tfault,deld,skew,prop,dmat,stn,histry,
     & s,stemp,gauss,ien,infin,lm,lmx,lmf,mat,iddmat,ngauss,nddmat,
     & ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,neq,numel,numnp,numfn,
     & numslp,numrot,nskdim,ipstrs,nstep,lgdefp,ibbarp,nhist,lastep,
     & idout,kto,kw,ivisc,iplas,imhist)
c
c...program to compute forces due to kinematic boundary conditions
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ngauss,nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,neq
      integer numel,numnp,numfn,numslp,numrot,nskdim,ipstrs,nstep,lgdefp
      integer ibbarp,nhist,lastep,idout,kto,kw,ivisc,iplas,imhist
      integer ien(nen,numel),infin(numel),lm(ndof,nen,numel)
      integer lmx(ndof,nen,numel),lmf(nen,numel),mat(numel)
      double precision b(neq),x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn),deld(ndof,numnp)
      double precision skew(nskdim,numnp),prop(nprop,numat)
      double precision dmat(nddmat,ngauss,ndmat),stn(nstr,ngauss,numel)
      double precision histry(nhist,lastep+1),s(nee*nee),stemp(nee*nee)
      double precision gauss(nsd+1,ngauss)
c
c...  included dimension and type statements
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
      integer ldtmp,io1,n,m,imat,i
      double precision p(24),dld(24),ptmp(30)
c
cdebug      write(6,*) "Hello from formdf_f!"
cdebug      write(6,*) "gauss:",((gauss(jdb,idb),jdb=1,nsd+1),idb=1,ngauss)
cdebug      write(6,*) "neq, b:",neq,(b(idb),idb=1,neq)
cdebug      write(6,*) "From formdf_f, ngauss, nsd, gauss: ",ngauss,nsd,
cdebug     & ((gauss(idb,jdb),idb=1,nsd+1),jdb=1,ngauss)
c
      ldtmp=lgdefp
      if(ipstrs.eq.ione.and.nstep.eq.izero) ldtmp=izero
      io1=izero
      if(ldtmp.gt.izero) io1=ione
      do n=1,numel
        m=mat(n)
        imat=n
        if((ivisc.eq.izero).and.(iplas.eq.izero)) imat=m
        call fill(p,zero,nee)
c
c...localize displacement boundary conditions
c
        call ldisbc(dld,deld,ien(1,n),lm(1,1,n),ndof,nen,numnp)
        do i=1,ndof*nen
          if(dld(i).ne.zero) go to 100
        end do
        go to 150
100     continue
c
c...  form element stiffness matrix
c
        call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &   idout,kto,kw,imhist)
        call formes(x,d,dx,tfault,dmat(1,1,imat),stn(1,1,n),skew,s,
     &   stemp,ptmp,gauss,ien(1,n),lmx(1,1,n),lmf(1,n),iddmat,infin(n),
     &   n,ngauss,nddmat,nprop,nsd,ndof,nstr,nen,nee,numnp,numfn,numslp,
     &   numrot,nskdim,io1,ibbarp,ldtmp,idout,kto,kw)
	call dsymv("u",nee,one,s,nee,dld,ione,zero,p,ione)
        call addfor(b,p,lm(1,1,n),lmx(1,1,n),neq,nee)
150     continue
      end do
cdebug      write(6,*) "from end of formdf, b:",(b(idb),idb=1,neq)
      return
      end
c
c version
c $Id: formdf.f,v 1.1 2004/06/15 20:03:35 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
