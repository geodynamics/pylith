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
      subroutine gload(
     & b,bres,x,d,dx,tfault,skew,grav,gvec1,gvec2,                      ! global
     & gauss,ien,lm,lmx,lmf,mat,infin,prop,histry,fulout,               ! elemnt
     & neq,nee,nsd,numnp,ndof,nen,ngauss,numel,numfn,numslp,nskdim,     ! dimens
     & numrot,nprop,numat,nhist,lastep,lgdefp,nstep,ipstrs,idebug,idout,! dimens
     & kto,kw,imhist)                                                   ! dimens
c
c...add body forces due to gravity
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer neq,nee,nsd,numnp,ndof,nen,ngauss,numel,numfn,numslp
      integer nskdim,numrot,nprop,numat,nhist,lastep,lgdefp,nstep,ipstrs
      integer idebug,idout,kto,kw,imhist
      integer ien(nen,numel),lm(ndof,nen,numel),lmx(ndof,nen,numel)
      integer lmf(nen,numel),mat(numel),infin(numel)
      double precision b(neq),bres(neq),x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp),tfault(ndof,numfn)
      double precision skew(nskdim,numnp),grav(ndof),gvec1(neq)
      double precision gvec2(neq),gauss(nsd+1,ngauss),prop(nprop,numat)
      double precision histry(nhist,lastep+1)
      logical fulout
c
c...  defined constants
c
      include "rconsts.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
cdebug      integer idb,jdb
      integer ldtmp,io2,i,npage,n,m
      double precision gsum,dif
      double precision p(24),xl(24),dl(24),ptmp(30)
c
c
c...clear body force array from last time step and check to see
c   whether body forces are being applied
c
cdebug      write(6,*) "Hello from gload_f!"
cdebug      write(6,*) "From gload_f, ngauss,nsd,gauss: ",ngauss,nsd,
cdebug     & ((gauss(idb,jdb),idb=1,nsd+1),jdb=1,ngauss)
c
      ldtmp=lgdefp
      if(ipstrs.eq.1.and.nstep.eq.0) ldtmp=0
      gsum=zero
      io2=0
      if(ldtmp.gt.0) io2=1
      do i=1,ndof
        gsum=gsum+grav(i)*grav(i)
      end do
      if(gsum.eq.zero) return
      call fill(gvec2,zero,neq)
      npage=50
c
c...loop over elements
c
      do n=1,numel
        m=mat(n)
        call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &   idout,kto,kw,imhist)
        call lcoord(x,xl,ien(1,n),nen,nsd,numnp)
        call fill(p,zero,nee)
c
c...update nodal positions and compute corresponding body force
c   store results in gvec2 (total body force vector)
c
        if(ldtmp.gt.0) call ldupdat(d,dx,tfault,dl,xl,ien(1,n),
     &   lmx(1,1,n),lmf(1,n),ndof,nsd,nen,numnp,numfn,numslp,io2,ldtmp)
        call gravldql(p,xl,grav,ptmp,gauss,ien(1,n),infin(n),n,nen,nee,
     &   nsd,ngauss,nprop,idout,kto,kw)
        if(numrot.ne.0) call rpforc(p,skew,ien(1,n),ndof,numnp,nen,
     &   nskdim)
        call addfor(gvec2,p,lm(1,1,n),lmx(1,1,n),neq,nee)
c
c...print out local load vectors if requested for debugging
c
        if(idebug.eq.1.and.idout.gt.1.and.fulout) then
          if(n.eq.1.or.mod(n,npage).eq.0) write(kw,1000)
          call prntforc(n,p,ien(1,n),nen,ndof,idout,kw)
        end if
      end do
c
c...find difference between computed total body force for new nodal
c   positions (gvec2) and current body force (gvec1) and update
c   body force vector and global load vector by this amount.
c
      do i=1,neq
        dif=gvec2(i)-gvec1(i)
        gvec1(i)=gvec2(i)
        b(i)=b(i)+dif
        bres(i)=bres(i)+dif
      end do
 1000 format(//' local forces computed by gload follow'//)
      return
      end
c
c version
c $Id: gload.f,v 1.1 2004/07/01 20:20:52 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
