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
      subroutine viscos(alnz,pcg,zcg,ja,                                ! sparse
     & b,btot,bres,pvec,gvec1,gvec2,                                    ! force
     & x,d,dx,deld,deldx,dprev,dcur,dxcur,id,idx,skew,histry,           ! global
     & ien,infin,mat,lm,lmx,lmf,prop,gauss,                             ! elemnt
     & ibond,bond,                                                      ! bc
     & dmat,stn,scur,st0,eps,deps,beta,dbeta,betb,dbetb,iddmat,         ! stress
     & ielno,iside,ihistry,pres,pdir,                                   ! tractn
     & maxstp,delt,alfa,maxit,maxitc,lgdef,ibbar,utol,ftol,etol,itmax,  ! timdat
     & iprint,                                                          ! output
     & fault,nfault,dfault,tfault,idftn,                                ! split
     & idslp,ipslp,diforc,idhist,                                       ! slip
     & iwink,wink,iwinkx,winkx,                                         ! wink
     & s,stemp,                                                         ! local
     & gcurr,gi,gprev,grav,gtol,ncodat,ndimens,                         ! info
     & npar,nprint,nsiter,nsysdat,ntimdat,nunits,nvisdat,rgiter,        ! info
     & rmin,rmult,rtimdat,                                              ! info
     & ofile,pfile)                                                     ! files

c
c...subroutine to solve the time dependent problem and perform the
c   time stepping
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ja(*),id(*),idx(*),ien(*),infin(*),mat(*),lm(*),lmx(*)
      integer lmf(*),ibond(*),ielno(*),iside(*),ihistry(*),maxstp(*)
      integer maxit(*),maxitc(*),lgdef(*),ibbar(*),itmax(*),iprint(*)
      integer nfault(*),idftn(*),idslp(*),ipslp(*),idhist(*),iwink(*)
      integer iwinkx(*)
      double precision alnz(*),pcg(*),zcg(*),b(*),btot(*),bres(*)
      double precision pvec(*),gvec1(*),gvec2(*),x(*),d(*),dx(*),deld(*)
      double precision deldx(*),dprev(*),dcur(*),dxcur(*),skew(*)
      double precision histry(*),prop(*),gauss(*),bond(*),dmat(*),stn(*)
      double precision scur(*),st0(*),eps(*),deps(*),beta(*),dbeta(*)
      double precision betb(*),dbetb(*),pres(*),pdir(*),delt(*),alfa(*)
      double precision utol(*),ftol(*),etol(*),fault(*),dfault(*)
      double precision tfault(*),diforc(*),wink(*),winkx(*),s(*)
      double precision stemp(*)
      character ofile*(*),pfile*(*)
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
      intrinsic abs,mod
c
c...  local variables
c
cdebug      integer idb
      integer indexx,ntot,jcyc,nfirst,naxstp,i,j,ii
      double precision time
      logical ltim,fulout,unlck,unlckf,skc,reform
c
c...  included variable definitions
c
      include "ndimens_def.inc"
      include "npar_def.inc"
      include "nprint_def.inc"
      include "nsysdat_def.inc"
      include "nunits_def.inc"
      include "nvisdat_def.inc"
      include "ntimdat_def.inc"
      include "rtimdat_def.inc"
c
c...  open output files for appending, if necessary
c
cdebug      write(6,*) "Hello from viscos_f!"
cdebug      write(6,*) "From viscos_f, ngauss,nsd,gauss:",ngauss,nsd,
cdebug     & (gauss(idb),idb=1,(nsd+1)*ngauss)
c
      if(idout.gt.1) open(kw,file=ofile,status="old",access="append")
      if(idsk.eq.0) open(kp,file=pfile,status="old",access="append")
      if(idsk.eq.1) open(kp,file=pfile,status="old",form="unformatted",
     & access="append")
c
c...signal user that viscous computation is begun
c
      write(kto,2000)
c
c...loop over complete cycles
c
      reform=.false.
      if(ireform.eq.1) reform=.true.
      indexx=1
      ntot=0
      do jcyc=1,ncycle
        if(ncycle.gt.1) write(kto,2001) jcyc
        nfirst=izero
        naxstp=izero
        nstep=izero
        time=zero
c
c...loop over time step groups
c
        do i=2,nintg
c
c...define constants to control stepping in current group
c
          call const(maxstp,delt,alfa,maxit,maxitc,lgdef,ibbar,utol,
     &     ftol,etol,itmax,nintg,i,naxstp,nfirst,rtimdat,deltp,alfap,
     &     ntimdat,nstep,maxitp,maxitcp,lgdefp,ibbarp,itmaxp,gtol)
cdebug          write(6,*) nintg,i,naxstp,nfirst,deltp,alfap,maxitp,maxitcp,
cdebug     &     lgdefp,ibbarp,itmaxp,(gtol(idb),idb=1,3)
          ltim=.true.
          if(alfap.eq.zero) ltim=.false.
c
c...loop over time steps in current group
c
          do j=nfirst,naxstp
            ntot=ntot+1
            nstep=nstep+1
            ntimdat(1)=nstep
            time=time+deltp
            skc=(numslp.ne.0.and.(iskopt.eq.2.or.(iskopt.le.0.and.
     &       abs(iskopt).eq.nstep)))
c
c...clear arrays at beginning of time step
c
            call fill(deld,zero,ndof*numnp)
            call fill(deldx,zero,ndof*numnp)
            if(skc) then
              call skclear(idslp,skew,numsn,nskdim,numnp)
              call skcomp(x,d,skew,idslp,ipslp,nsd,ndof,nskdim,npdim,
     &         ipstrs,numsn,numnp,nstep,lgdefp,kto)
            end if
c
c...test for reform and refactor interval, whether full output
c   occurs in this step
c
            if(nwink.ne.0) call cklock(iwink,histry,ltim,nwink,nstep,
     &       nhist,lastep,unlck)
            if(nwinkx.ne.0) call cklock(iwinkx,histry,ltim,nwinkx,nstep,
     &       nhist,lastep,unlckf)
            reform=(alfap.ne.0.0.or.lgdefp.ne.0).and.
     &       (mod(j,maxitp).eq.0).or.ltim
            ireform=izero
            if(reform) ireform=1
            ntimdat(10)=ireform
            fulout=.false.
            if(ntot.eq.iprint(indexx)) fulout=.true.
c
            if(idout.gt.1) write(kw,1000) time,ntot,jcyc
C***********************************
            if(fulout.and.idsk.eq.0) write(kp,700) ntot
C***********************************
            if(fulout.and.idsk.eq.0) write(kp,4000) time
            if(fulout.and.idsk.eq.1) write(kp) ntot
            if(fulout.and.idsk.eq.1) write(kp) time
            write(kto,5000) time,ntot,lastep*ncycle
c*            call flush(kto)
c*            call flush(kw)
c*            call flush(kp)
c
c...apply boundary conditions
c
            call load(id,ibond,bond,d,deld,btot,histry,deltp,ndof,numnp,
     &       neq,nhist,nstep,lastep,idout,kto,kw)
c
c...compute current split node displacements
c
            if(numfn.ne.0) call loadf(fault,dfault,histry,deltp,nfault,
     &       nstep,ndof,numfn,nhist,lastep,idout,kto,kw)
c
c...add loads from changes in differential forces across internal
c        interfaces
c
            if(numdif.ne.0) call loadx(btot,diforc,histry,idx,idhist,
     &       ndof,neq,numnp,nhist,nstep,lastep,idout,kto,kw)
c
c...compute change in load vector if winkler forces are removed
c
            call fill(zcg,zero,neq)
            if(nwink.ne.0.and.unlck) call unlock(zcg,btot,id,iwink,
     &       idhist,ibond,bond,histry,nstep,ndof,numnp,nwink,nhist,neq,
     &       numdif,lastep,ione)
            if(nwinkx.ne.0.and.unlckf) call unlock(zcg,btot,idx,iwinkx,
     &       idhist,ibond,diforc,histry,nstep,ndof,numnp,nwinkx,nhist,
     &       neq,numdif,lastep,itwo)
	    call daxpy(neq,one,zcg,ione,btot,ione)
c
c...compute forces due to applied displacements and split nodes
c
            ii=ione
            call formmat(stn,dmat,deps,beta,betb,
     &       prop,mat,iddmat,
     &       histry,rtimdat,ntimdat,rgiter,ii,
     &       ngauss,nddmat,ndmat,nprop,numat,ndof,nstr,numel,ipstrs,
     &       nhist,lastep,idout,kto,kw)
            call formdf(b,x,d,dx,tfault,deld,skew,prop,dmat,stn,histry,
     &       s,stemp,gauss,ien,infin,lm,lmx,lmf,mat,iddmat,ngauss,
     &       nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,neq,numel,
     &       numnp,numfn,numslp,numrot,nskdim,ipstrs,nstep,lgdefp,
     &       ibbarp,nhist,lastep,idout,kto,kw,ivisc,iplas,imhist)
            if(numfn.ne.0) call formf(b,x,d,dx,skew,histry,
     &       ien,lm,lmx,lmf,gauss,
     &       mat,infin,prop,dmat,stn,
     &       nfault,dfault,tfault,
     &       s,stemp,iddmat,
     &       ngauss,nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,neq,
     &       numel,numnp,numfn,numslp,numrot,nskdim,ipstrs,nstep,lgdefp,
     &       ibbarp,nhist,lastep,idout,kto,kw,ivisc,iplas,imhist)
c
c...compute iterative solution to the elastic problem
c   return if convergence achieved or maximum iterations exceeded
c
            call iterate(alnz,pcg,zcg,ja,                               ! sparse
     &       b,btot,bres,pvec,gvec1,gvec2,                              ! force
     &       x,d,dx,deld,deldx,dprev,dcur,dxcur,id,idx,skew,histry,     ! global
     &       ien,infin,mat,lm,lmx,lmf,prop,gauss,                       ! elemnt
     &       dmat,stn,scur,st0,eps,deps,beta,dbeta,betb,dbetb,iddmat,   ! stress
     &       ielno,iside,ihistry,pres,pdir,                             ! tractn
     &       nfault,dfault,tfault,                                      ! split
     &       idslp,ipslp,                                               ! slip
     &       iwink,wink,iwinkx,winkx,                                   ! wink
     &       s,stemp,                                                   ! local
     &       gcurr,gi,gprev,grav,gtol,ncodat,                           ! info
     &       ndimens,npar,nprint,nsiter,nsysdat,ntimdat,nunits,nvisdat, ! info
     &       rgiter,rmin,rmult,rtimdat)                                 ! info
c
c...print displacements at all nodes when requested.
c
            if(fulout) then
              call printd(d,deld,deltp,idslp,ndof,numnp,numnp,ione,
     &         idout,idsk,kto,kw,kp)
              call printf(tfault,dfault,deltp,nfault,ndof,numfn,idout,
     &         idsk,kw,kp)
              call printd(dx,deldx,deltp,idslp,ndof,numnp,numsn,itwo,
     &         idout,idsk,kto,kw,kp)
              call printl(idx,iwinkx,idslp,histry,ndof,numsn,numnp,
     &         nstep,nhist,nwinkx,lastep,idsk,kp)
            end if
c
c...print stresses in all elements when requested
c
            if(fulout) call prints(stn,eps,deps,beta,dbeta,betb,dbetb,
     &       nstr,ngauss,numel,nstep,idout,idsk,kw,kp,ivisc,iplas)
            ltim=.false.
            if(fulout) indexx=indexx+1
            if(indexx.gt.icontr) indexx=icontr
          end do
        end do
      end do
      write(kto,800) ntimdat(7),ntimdat(8),ntimdat(9)
      if(idout.gt.1) write(kw,800) ntimdat(7),ntimdat(8),ntimdat(9)
      if(idout.gt.1) close(kw)
      close(kp)
c
 700  format('STEP ',i7)
 800  format(/," Total number of equilibrium iterations        = ",i7,/,
     &         " Total number of stiffness matrix reformations = ",i7,/,
     &         " Total number of displacement subiterations    = ",i7)
 1000 format(////' output which follows is at time= ',1pe15.8/
     & ' step # ',i7,' in cycle # ',i7//)
 2000 format(/'Time-dependent solution is begun:'/)
 2001 format('     working on cycle ',i7)
 4000 format(1pe15.8)
 5000 format(///,'   Working on time ',1pe15.8,', timestep #',i7,
     & '  out of',i7)
      return
      end
c
c version
c $Id: viscos.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
