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
      subroutine addpr(
     & b,bres,x,d,dx,tfault,histry,skew,                                ! global
     & ien,infin,lm,lmx,lmf,                                            ! elemnt
     & ielno,iside,ihstry,pres,pdir,pvec,gvec2,fulout,                  ! press
     & nen,numnp,neq,nee,numrot,lastep,nhist,                           ! dimens
     & nstep,lgdefp,numel,numpr,numfn,numslp,ipstrs,idout,idebug,kto,kw)! dimens
c
c...subroutine to add pressures to load vector
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
c
c...  subroutine arguments
c
      integer nen,numnp,neq,nee,numrot,lastep
      integer nhist,nstep,lgdefp,numel,numpr,numfn,numslp,ipstrs,idout
      integer idebug,kto,kw
      integer ien(nen,numel),infin(numel),lm(ndof,nen,numel)
      integer lmx(ndof,nen,numel),lmf(nen,numel),ielno(numpr)
      integer iside(numpr),ihstry(numpr)
      double precision b(neq),bres(neq),x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp),tfault(ndof,numfn)
      double precision histry(nhist,lastep+1),skew(nskdim,numnp)
      double precision pres(nen/2,numpr),pdir(npdir,numpr),pvec(neq)
      double precision gvec2(neq)
      logical fulout
c
c...  defined constants
c
      include "rconsts.inc"
c
c...  local variables
c
      double precision diff,dif
      double precision dl(24),p(24),xl(24),pload(4)
      integer io2,npage,ldtmp,k,ihist,n,i
c
c...clear gvec2 so that pressure contributions to load vector can be
c   stored there for constant pressure BC
c
cdebug      write(6,*) "Hello from addpr_f!"
c
      call fill(gvec2,zero,neq)
      io2=1
      npage=50
      ldtmp=lgdefp
      if(ipstrs.eq.1.and.nstep.eq.0) ldtmp=0
c
c...loop over pressure bc
c
      do k=1,numpr
        ihist=ihstry(k)
        if(ihist.gt.nhist) then
          if(idout.gt.1) write(kw,*)
     &     ' attempt to use undefined load history # ',ihist
          write(kto,*) ' attempt to use undefined load history # ',ihist
          stop
        end if
        if(ihist.ne.0.or.nstep.eq.0.or.lgdefp.ne.0) then
          if(ihist.eq.0) then
            pload(1)=pres(1,k)
            pload(2)=pres(2,k)
            if(nsd.eq.3) then
              pload(3)=pres(3,k)
              pload(4)=pres(4,k)
            end if
          else if(ihist.gt.0) then
            diff=histry(ihist,nstep+1)
            if(nstep.gt.0.and.lgdefp.eq.0) diff=diff-histry(ihist,nstep)
            pload(1)=pres(1,k)*diff
            pload(2)=pres(2,k)*diff
            if(nsd.eq.3) then
              pload(3)=pres(3,k)*diff
              pload(4)=pres(4,k)*diff
            end if
          end if
          n=ielno(k)
          call lcoord(x,xl,ien(1,n),nen,numnp)
c
c...update nodal positions for large deformation formalism
c
          if(ldtmp.ge.1) call ldupdat(d,dx,tfault,dl,xl,ien(1,n),
     &     lmx(1,1,n),lmf(1,n),nen,numnp,numfn,numslp)
c
c...compute local load vector and add it to global vector or
c   temporary vector if constant pressure bc are desired for large
c   deformations
c
          call presurql(pload,pdir(1,k),xl,p,ien(1,n),iside(k),infin(n),
     &     n,nen,idout,kto,kw)
          if(numrot.ne.0) call rpforc(p,skew,ien(1,n),numnp,nen)
          if(lgdefp.eq.0) call addfor(b,p,lm(1,1,n),lmx(1,1,n),neq,nee)
          if(lgdefp.eq.0) call addfor(bres,p,lm(1,1,n),lmx(1,1,n),neq,
     &     nee)
          if(lgdefp.gt.0) call addfor(gvec2,p,lm(1,1,n),lmx(1,1,n),neq,
     &     nee)
c
c...print out local load vectors if requested for debugging
c
          if(idebug.eq.1.and.idout.gt.1.and.fulout) then
            if(n.eq.1.or.mod(n,npage).eq.0) write(kw,1000)
            call prntforc(n,p,ien(1,n),nen,idout,kw)
          end if
        end if
      end do
c
c...find difference between pressure loads from last time step and
c   loads for current geometry for constant pressure BC.
c
      if(lgdefp.gt.0) then
        do i=1,neq
          dif=gvec2(i)-pvec(i)
          pvec(i)=gvec2(i)
          b(i)=b(i)+dif
          bres(i)=bres(i)+dif
        end do
      end if
 1000 format(//' local forces computed by addpr follow'//)
      return
      end
c
c version
c $Id: addpr.f,v 1.2 2004/08/12 01:03:16 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
