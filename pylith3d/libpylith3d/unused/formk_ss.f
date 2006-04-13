c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
c  All worldwide rights reserved.  A license to use, copy, modify and
c  distribute this software for non-commercial research purposes only
c  is hereby granted, provided that this copyright notice and
c  accompanying disclaimer is not modified or removed from the software.
c
c  DISCLAIMER:  The software is distributed "AS IS" without any express
c  or implied warranty, including but not limited to, any implied
c  warranties of merchantability or fitness for a particular purpose
c  or any warranty of non-infringement of any current or pending patent
c  rights.  The authors of the software make no representations about
c  the suitability of this software for any particular purpose.  The
c  entire risk as to the quality and performance of the software is with
c  the user.  Should the software prove defective, the user assumes the
c  cost of all necessary servicing, repair or correction.  In
c  particular, neither Rensselaer Polytechnic Institute, nor the authors
c  of the software are liable for any indirect, special, consequential,
c  or incidental damages related to the software, to the maximum extent
c  the law permits.
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
