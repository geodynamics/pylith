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
c $Id: formdf.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
