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
c
      subroutine formk_ss(
     & alnz,ja,nnz,neq,                                                 ! sparse
     & id,x,d,iwink,wink,                                               ! global
     & tfault,numfn,
     & idx,dx,iwinkx,winkx,
     & skew,numrot,
     & histry,nhist,nstep,lastep,
     & ien,lm,lmx,lmf,mat,infin,prop,gauss,                             ! elemnt
     & dmat,stn,eps,beta,betb,iter,                                     ! stress
     & s,stemp,iddmat,                                                  ! local
     & rtimdat,ntimdat,rgiter,                                          ! timdat
     & ngauss,nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,! dimens
     & numel,numnp,numfn,numslp,numsn,nskdim,ipstrs,             ! dimens
     & nwink,nwinkx,idout,kto,kw)                          ! dimens
c
c...program to form the stiffness matrix k
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iter,ngauss,nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee
      integer nnz,neq,numel,numnp,numfn,numslp,numsn,numrot,nskdim
      integer ipstrs,nwink,nwinkx,nhist,lastep,idout,kto,kw
      integer ja(nnz),id(ndof,numnp),idx(ndof,numnp),iwink(2,nwink)
      integer iwinkx(2,nwinkx),ien(nen,numel),lm(ndof,nen,numel)
      integer lmx(ndof,nen,numel),lmf(nen,numel),mat(numel),infin(numel)
      double precision alnz(nnz),x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp),tfault(ndof,numfn)
      double precision skew(nskdim,numnp),wink(nwink),winkx(nwinkx)
      double precision histry(nhist,lastep+1),prop(nprop,numat)
      double precision gauss(nsd+1,ngauss),dmat(nddmat,ngauss,ndmat)
      double precision stn(nstr,ngauss,numel),eps(nstr,ngauss,numel)
      double precision beta(nstr,ngauss,numel),betb(nstr,ngauss,numel)
      double precision s(nee*nee),stemp(nee*nee)
c
c...  included dimension and type statements
c
      include "iddmat_dim.inc"
      include "ntimdat_dim.inc"
      include "rgiter_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  local variables
c
cdebug      integer idb,jdb
      integer ldtmp,io1,m,l,n,imat,nzero
      double precision ptmp(30),dummy(1)
      logical*4 newmat
c
c...  included variable definitions
c
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from formk_f!"
cdebug      write(6,*) "From formk_f, ngauss,nsd,gauss: ",ngauss,nsd,
cdebug     & ((gauss(jdb,idb),jdb=1,nsd+1),idb=1,ngauss)
c
      call fill(alnz,zero,nnz)
      nrftot=nrftot+1
      ntimdat(8)=nrftot
      ldtmp=lgdefp
      if(ipstrs.eq.ione.and.nstep.eq.izero) ldtmp=izero
      io1=izero
      if(ldtmp.gt.izero) io1=ione
      newmat=.false.
      if(((ivisc.eq.ione).or.(iplas.eq.ione)).and.
     & nstep.gt.izero.and.iter.gt.ione) newmat=.true.
cdebug      write(6,*) ivisc,iplas,iter,newmat
      if((ivisc.eq.izero).and.(iplas.eq.izero)) then
        do m=1,numat
          call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &     idout,kto,kw,imhist)
          do l=1,ngauss
            call matinit(stn(1,l,m),eps(1,l,m),beta(1,l,m),betb(1,l,m),
     &       dmat(1,l,m),ptmp,rtimdat,iddmat,nstr,nddmat,nprop,
     &       ndof,ipstrs,nstep,lgdefp,ivisc,iplas)
          end do
        end do
      end if
c
c...loop over material groups
c
      do n=1,numel
        m=mat(n)
        call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &     idout,kto,kw,imhist)
        imat=n
        if((ivisc.eq.izero).and.(iplas.eq.izero)) imat=m
c
c...construct the local material matrix and insert it in the global
c   array dmat(nddmat,ngauss,ndmat)
c
        if(newmat) then
          do l=1,ngauss
            call matprtb(stn(1,l,n),eps(1,l,n),beta(1,l,n),betb(1,l,n),
     &       dmat(1,l,n),ptmp,rtimdat,rgiter,iddmat,n,nstr,ndof,
     &       nprop,nddmat,ipstrs,nstep,lgdefp,idout,kto,kw,ivisc,iplas)
cdebug            if(n.eq.20) write(6,*) (dmat(jdb,l,imat),jdb=1,nddmat)
          end do
        end if
c
c...construct the element stiffness matrix
c
        call formes(x,d,dx,tfault,dmat(1,1,imat),stn(1,1,n),skew,s,
     &   stemp,ptmp,gauss,ien(1,n),lmx(1,1,n),lmf(1,n),iddmat,infin(n),
     &   n,ngauss,nddmat,nprop,nsd,ndof,nstr,nen,nee,numnp,numfn,numslp,
     &   numrot,nskdim,io1,ibbarp,ldtmp,idout,kto,kw)
cdebug        if(n.eq.20) write(6,*) (s(jdb),jdb=1,100)
        call addstf(alnz,s,lm(1,1,n),lmx(1,1,n),ja,nee,
     &   numel,numsn,nnz)
      end do
c
c...add Winkler elements to stiffness matrix diagonals
c
      if(nwink.ne.0) call winklr(alnz,iwink,wink,histry,nstep,nwink,
     & nhist,nnz,lastep,idout,kto,kw)
      if(nwinkx.ne.0) call winklr(alnz,iwinkx,winkx,histry,nstep,nwinkx,
     & nhist,nnz,lastep,idout,kto,kw)
c
c...check stiffness matrix for zero or negative diagonals, and stop if
c   they are found, after printing out a list of the diagonals
c
      call ckdiag(alnz,nzero,neq,nnz,idout,kto,kw)
      if(nzero .eq. ione) then
        call printv(alnz,dummy,id,idx,neq,ndof,numnp,itwo,idout,kw)
        stop
      end if
cdebug      write(6,*) (alnz(jdb),jdb=1,neq)
      return
      end
c
c version
c $Id: formk_ss.f,v 1.1 2004/07/02 18:42:02 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
