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
      subroutine iterate(alnz,pcg,zcg,ja,                               ! sparse
     & b,btot,bres,pvec,gvec1,gvec2,                                    ! force
     & x,d,dx,deld,deldx,dprev,dcur,dxcur,id,idx,skew,histry,           ! global
     & ien,infin,mat,lm,lmx,lmf,prop,gauss,                             ! elemnt
     & dmat,stn,scur,st0,eps,deps,beta,dbeta,betb,dbetb,iddmat,         ! stress
     & ielno,iside,ihistry,pres,pdir,                                   ! tractn
     & nfault,dfault,tfault,                                            ! split
     & idslp,ipslp,                                                     ! slip
     & iwink,wink,iwinkx,winkx,                                         ! wink
     & s,stemp,                                                         ! local
     & gcurr,gi,gprev,grav,gtol,ncodat,ndimens,                         ! info
     & npar,nprint,nsiter,nsysdat,ntimdat,nunits,nvisdat,rgiter,        ! info
     & rmin,rmult,rtimdat)                                              ! info
c
c...subroutine to loop over iterations within each time step
c
c       This is a program segment of tecton that performs an iterative
c       solution until both the force tolerance and the energy tolerance
c       are achieved, or until the maximum number of iterations has
c       been exceeded.  In the case of large deformations, pressure BC
c       and body forces are re-evaluated over the new geometry.
c       This routine is called by both the elastic and viscoelastic
c       portions of the code.  When called from the elastic portion,
c       the elastic material matrix is used to form the stiffness matrix
c       and to evaluate the stresses.  When called from the viscoelastic
c       portion of the code, the effective stress function algorithm
c       is used.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ja(*),id(*),idx(*),ien(*),infin(*),mat(*),lm(*),lmx(*)
      integer lmf(*),ielno(*),iside(*),ihistry(*),nfault(*),idslp(*)
      integer ipslp(*),iwink(*),iwinkx(*)
      double precision alnz(*),pcg(*),zcg(*),b(*),btot(*),bres(*)
      double precision pvec(*),gvec1(*),gvec2(*),x(*),d(*),dx(*),deld(*)
      double precision deldx(*),dprev(*),dcur(*),dxcur(*),skew(*)
      double precision histry(*),prop(*),gauss(*),dmat(*),stn(*),scur(*)
      double precision st0(*),eps(*),deps(*),beta(*),dbeta(*),betb(*)
      double precision dbetb(*),pres(*),pdir(*),dfault(*),tfault(*)
      double precision wink(*),winkx(*),s(*),stemp(*)
c
c...  included dimension and type statements
c
      include "iddmat_dim.inc"
      include "gcurr_dim.inc"
      include "gi_dim.inc"
      include "gprev_dim.inc"
      include "grav_dim.inc"
      include "gtol_dim.inc"
      include "ncodat_dim.inc"
      include "ndimens_dim.inc"
      include "npar_dim.inc"
      include "nprint_dim.inc"
      include "nsiter_dim.inc"
      include "nsysdat_dim.inc"
      include "ntimdat_dim.inc"
      include "nunits_dim.inc"
      include "nvisdat_dim.inc"
      include "rgiter_dim.inc"
      include "rmin_dim.inc"
      include "rmult_dim.inc"
      include "rtimdat_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  intrinsic functions
c
      intrinsic mod,abs
c
c...  local variables
c
cdebug      integer idb
      integer i
      double precision dummy(1)
      logical fulout,converge,updats,skc,used,reform
c
c...  included variable definitions
c
      include "ncodat_def.inc"
      include "ndimens_def.inc"
      include "npar_def.inc"
      include "nprint_def.inc"
      include "nsiter_def.inc"
      include "nsysdat_def.inc"
      include "ntimdat_def.inc"
      include "nunits_def.inc"
      include "nvisdat_def.inc"
      include "rgiter_def.inc"
c
c...initialize convergence criteria
c
cdebug      write(6,*) "Hello from iterate_f!"
cdebug      write(6,*) "From iterate, ngauss,nsd,gauss: ",ngauss,nsd,
cdebug     & (gauss(idb),idb=1,(nsd+1)*ngauss)
c
      fulout=.true.
      converge=.false.
      reform=.false.
      if(ireform.eq.1) reform=.true.
      skc=.false.
      do i=1,3
        gcurr(i)=10.0d0*gtol(i)
        gprev(i)=10.0d0*gtol(i)
      end do
c
c...loop over iterations
c
      do i=1,itmaxp
        updats=((ipstrs.ne.1).or.(ipstrs.eq.1.and.nstep.eq.0.and.
     &   i.eq.1).or.(nstep.gt.0)).and.lgdefp.ge.1
        skc=i.gt.1.and.lgdefp.ge.1.and.(numslp.ne.0.and.(iskopt.eq.2.or.
     &   (iskopt.le.0.and.abs(iskopt).eq.nstep)))
        nittot=nittot+1
        ntimdat(7)=nittot
        reform=reform.or.(mod(i,maxitcp).eq.0)
        ireform=0
        if(reform) ireform=1
        ntimdat(10)=ireform
        used=nstep.gt.0.and.(nsol.eq.2.or.nsol.eq.4).and.i.eq.1
        if(i.gt.1) fulout=.false.
c
c...add pressure forces,if present, to global load vector
c
        if(numpr.ne.0.and.(i.eq.1.or.updats)) call addpr(
     &   btot,bres,x,d,dx,tfault,histry,skew,
     &   ien,infin,lm,lmx,lmf,
     &   ielno,iside,ihistry,pres,pdir,pvec,gvec2,fulout,
     &   nsd,ndof,nen,nskdim,npdir,numnp,neq,nee,numrot,lastep,nhist,
     &   nstep,lgdefp,numel,numpr,numfn,numslp,ipstrs,idout,idebug,kto,
     &   kw)
c
c...add gravity body forces to global load vector
c
        if(i.eq.1.or.updats) call gload(
     &   btot,bres,x,d,dx,tfault,skew,grav,gvec1,gvec2,
     &   gauss,ien,lm,lmx,lmf,mat,infin,prop,histry,fulout,
     &   neq,nee,nsd,numnp,ndof,nen,ngauss,numel,numfn,numslp,nskdim,
     &   numrot,nprop,numat,nhist,lastep,lgdefp,nstep,ipstrs,idebug,
     &   idout,kto,kw,imhist)
c
c...reform the stiffness matrix, if required
c
        if(reform) then
          if(skc) call skcomp(x,d,skew,idslp,ipslp,nsd,ndof,nskdim,
     &     npdim,ipstrs,numsn,numnp,nstep,lgdefp,kto)
          call formk(alnz,ja,
     &     id,idx,x,d,dx,tfault,skew,iwink,wink,iwinkx,winkx,histry,
     &     ien,lm,lmx,lmf,mat,infin,prop,gauss,
     &     dmat,stn,deps,beta,betb,i,
     &     s,stemp,iddmat,
     &     rtimdat,ntimdat,rgiter,
     &     ngauss,nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,nnz,
     &     neq,numel,numnp,numfn,numslp,numsn,numrot,nskdim,ipstrs,
     &     nwink,nwinkx,nhist,lastep,idout,kto,kw)
        end if
c
c...if icode .eq. 1 print the stiffness matrix diagonals and stop here
c
        if(icode.eq.1.and.nstep.eq.0.and.i.eq.1) then
          call printv(alnz,dummy,id,idx,neq,ndof,numnp,izero,idout,kw)
          stop
        else if(icode.eq.1) then
          stop
        end if
c
c...for first iteration compute residual force vector
c
        if(i.eq.1) call bdiff(b,btot,bres,neq)
c
c...compute the global displacement increment vector using a
c   preconditioned conjugate gradients iterative solver.  Upon
c   return the vector gvec2 contains the displacements.
c
        call pcginv(alnz,bres,gvec2,b,pcg,zcg,dprev,rmin,rmult,
     &   gcurr,gprev,ja,nsiter,neq,nnz,ndtot,idout,kto,kw,used)
        ntimdat(9)=ndtot
        if(nsol.eq.2.or.nsol.eq.4) then
          if(i.eq.1) call fill(dprev,zero,neq)
          call daxpy(neq,one,gvec2,ione,dprev,ione)
        end if
c
c...for first iteration, update displacements to reflect boundary
c   conditions
c
        if(i.eq.1) then
          if(numfn.ne.0.and.numrot.ne.0) call rsplit(nfault,dfault,
     &     skew,ndof,numfn,numnp,nskdim)
          if(numfn.ne.0) call daxpy(ndof*numfn,one,dfault,ione,tfault,
     &     ione)
          if(nstep.eq.0) then
            if(numrot.ne.0) call rdisp(dcur,skew,ndof,numnp,nskdim)
            call dcopy(ndof*numnp,dcur,ione,d,ione)
            call dcopy(ndof*numnp,dcur,ione,deld,ione)
          else
            if(numrot.ne.0) call rdisp(deld,skew,ndof,numnp,nskdim)
            call daxpy(ndof*numnp,one,deld,ione,d,ione)
          end if
        end if
c
c...localize the displacement increments in dcur(ndof,numnp) and
c   dxcur(ndof,numnp)
c
        call fill(dcur,zero,ndof*numnp)
        call fill(dxcur,zero,ndof*numnp)
        call disp(gvec2,dcur,id,ndof,numnp,neq)
        if(numslp.ne.0) call disp(gvec2,dxcur,idx,ndof,numnp,neq)
c
c...rotate skewed coordinates to global system
c
        if(numrot.ne.0) then
          call rdisp(dcur,skew,ndof,numnp,nskdim)
          if(numslp.ne.0) call rdisp(dxcur,skew,ndof,numnp,nskdim)
        end if
c
c...compute contribution to global load vector from winkler boundary
c   conditions
c
        if(nwink.ne.0) call winklf(btot,gvec2,iwink,wink,histry,nwink,
     &   nhist,nstep,neq,lastep)
        if(nwinkx.ne.0) call winklf(btot,gvec2,iwinkx,winkx,histry,
     &   nwinkx,nhist,nstep,neq,lastep)
c
c...update the total displacement and the displacement increment
c   for this step
c
        call daxpy(ndof*numnp,one,dcur,ione,d,ione)
        call daxpy(ndof*numnp,one,dcur,ione,deld,ione)
        if(numslp.ne.0) call daxpy(ndof*numnp,one,dxcur,ione,dx,ione)
        if(numslp.ne.0) call daxpy(ndof*numnp,one,dxcur,ione,deldx,ione)
c
c...integrate the stresses and compute the equivalent nodal loads
c
cdebug        write(6,*) "From iterate right before stresn call:"
cdebug        write(6,*) "gauss:",(gauss(idb),idb=1,(nsd+1)*ngauss)
        call stresn(x,b,d,dx,tfault,stn,deps,beta,betb,scur,st0,dbeta,
     &   dbetb,skew,ien,lm,lmx,lmf,dmat,mat,prop,histry,infin,gauss,
     &   rtimdat,stol,iddmat,nen,numel,ndof,nsd,numnp,neq,nee,
     &   nstr,ngauss,nppts,ngem,nskdim,nhist,nprop,numat,numfn,numslp,
     &   numrot,lastep,nstep,lgdefp,ibbarp,ivisc,iplas,imhist,ipstrs,
     &   nprestr,nddmat,ndmat,idebug,idout,kto,kw,fulout)
c
c...compute the out-of-balance forces and convergence criteria
c
        call residu(b,bres,btot,gvec2,gtol,gi,gprev,gcurr,id,idx,neq,
     &   ndof,numnp,i,itmaxp,idebug,idout,kto,kw,converge)
c
c...if solution has converged, set equilibrium stresses and creep
c   strains to their current values
c
        if(converge) then
          call dcopy(nstr*numel*ngauss,scur,ione,stn,ione)
          if(ivisc.eq.1) call daxpy(nstr*numel*ngauss,one,dbeta,ione,
     &     beta,ione)
          if(iplas.eq.1) call daxpy(nstr*numel*ngauss,one,dbetb,ione,
     &     betb,ione)
          call daxpy(nstr*numel*ngauss,-one,eps,ione,deps,ione)
          call daxpy(nstr*numel*ngauss,one,deps,ione,eps,ione)
          return
        end if
        reform=.false.
c*        call flush(kto)
c*        call flush(kw)
c*        call flush(kp)
      end do
      return
      end
c
c version
c $Id: iterate.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
