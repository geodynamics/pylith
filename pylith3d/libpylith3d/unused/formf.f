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
      subroutine formf(
     & b,x,d,dx,skew,histry,                                             ! global
     & ien,lm,lmx,lmf,gauss,                                             ! elemnt
     & mat,infin,prop,dmat,stn,                                          ! stress
     & nfault,dfault,tfault,                                             ! split
     & s,stemp,iddmat,                                                   ! local
     & ngauss,nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,neq,numel,  ! dimens
     & numnp,numfn,numslp,numrot,nskdim,ipstrs,nstep,lgdefp,ibbarp,      ! dimens
     & nhist,lastep,idout,kto,kw,ivisc,iplas,imhist)                     ! dimens
c
c      generates forces due to faulted nodes
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ngauss,nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,neq
      integer numel,numnp,numfn,numslp,numrot,nskdim,ipstrs,nstep,lgdefp
      integer ibbarp,nhist,lastep,idout,kto,kw,ivisc,iplas,imhist
      integer ien(nen,numel),lm(ndof,nen,numel),lmx(ndof,nen,numel)
      integer lmf(nen,numel),mat(numel),infin(numel),nfault(3,numfn)
      double precision b(neq),x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision skew(nskdim,numnp),histry(nhist,lastep+1)
      double precision gauss(nsd+1,ngauss),prop(nprop,numat)
      double precision dmat(nddmat,ngauss,ndmat),stn(nstr,ngauss,numel)
      double precision dfault(ndof,numfn),tfault(ndof,numfn),s(nee*nee)
      double precision stemp(nee*nee)
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
      integer ldtmp,io1,i,n,m,imat
      double precision p(24),dl(24),ptmp(30)
c
cdebug      write(6,*) "Hello from formf_f!"
cdebug      write(6,*) "neq, b:",neq,(b(idb),idb=1,neq)
cdebug      write(6,*) "From formf_f, ngauss,nsd,gauss: ",ngauss,nsd,
cdebug     & ((gauss(jdb,idb),jdb=1,nsd+1),idb=1,ngauss)
c
      ldtmp=lgdefp
      if(ipstrs.eq.ione.and.nstep.eq.izero) ldtmp=izero
      io1=izero
      if(ldtmp.gt.izero) io1=ione
      do i = 1,numfn
        n=nfault(1,i)
        m=mat(n)
        imat=n
        if((ivisc.eq.izero).and.(iplas.eq.izero)) imat=m
        call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &   idout,kto,kw,imhist)
        call formes(x,d,dx,tfault,dmat(1,1,imat),stn(1,1,n),skew,s,
     &   stemp,ptmp,gauss,ien(1,n),lmx(1,1,n),lmf(1,n),iddmat,infin(n),
     &    n,ngauss,nddmat,nprop,nsd,ndof,nstr,nen,nee,numnp,numfn,
     &    numslp,numrot,nskdim,io1,ibbarp,ldtmp,idout,kto,kw)
        call lflteq(dl,dfault(1,i),nfault(1,i),ien(1,n),nen,ndof)
	call dsymv("u",nee,one,s,nee,dl,ione,zero,p,ione)
        call addfor(b,p,lm(1,1,n),lmx(1,1,n),neq,nee)
      end do
cdebug      write(6,*) "form end of formf, b:",(b(idb),idb=1,neq)
      return
      end
c
c version
c $Id: formf.f,v 1.1 2004/06/21 20:59:05 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
