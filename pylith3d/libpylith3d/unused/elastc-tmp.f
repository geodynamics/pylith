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
      subroutine elastc(alnz,pcg,zcg,ja,                                ! sparse
     & b,btot,bres,pvec,gvec1,gvec2,                                    ! force
     & x,d,dx,deld,deldx,dprev,dcur,dxcur,id,idx,skew,histry,           ! global
     & ien,infin,mat,lm,lmx,lmf,prop,gauss,                             ! elemnt
     & ibond,bond,                                                      ! bc
     & dmat,stn,scur,st0,eps,deps,beta,dbeta,betb,dbetb,iddmat,         ! stress
     & ielno,iside,ihistry,pres,pdir,                                   ! tractn
     & maxstp,delt,alfa,maxit,maxitc,lgdef,ibbar,utol,ftol,etol,itmax,  ! timdat
     & fault,nfault,dfault,tfault,idftn,                                ! split
     & idslp,ipslp,diforc,idhist,                                       ! slip
     & iwink,wink,iwinkx,winkx,                                         ! wink
     & s,stemp,                                                         ! local
     & gcurr,gi,gprev,grav,gtol,ncodat,ndimens,                         ! info
     & npar,nprint,nsiter,nsysdat,ntimdat,nunits,nvisdat,rgiter,        ! info
     & rmin,rmult,rtimdat,                                              ! info
     & ofile,pfile)                                                     ! files
c
c...subroutine to construct and solve the elastic problem
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer ja(*),id(*),idx(*),ien(*),infin(*),mat(*),lm(*),lmx(*)
      integer lmf(*),ibond(*),ielno(*),iside(*),ihistry(*),maxstp(*)
      integer maxit(*),maxitc(*),lgdef(*),ibbar(*),itmax(*),nfault(*)
      integer idftn(*),idslp(*),ipslp(*),idhist(*),iwink(*),iwinkx(*)
      double precision alnz(*),pcg(*),zcg(*),b(*),btot(*),bres(*),pvec(*),gvec1(*)
      double precision gvec2(*),x(*),d(*),dx(*),deld(*),deldx(*),dprev(*),dcur(*)
      double precision dxcur(*),skew(*),histry(*)
      double precision prop(*),gauss(*),bond(*),dmat(*),stn(*),scur(*),st0(*)
      double precision eps(*),deps(*),beta(*),dbeta(*),betb(*),dbetb(*),pres(*)
      double precision pdir(*),delt(*),alfa(*),utol(*),ftol(*),etol(*),fault(*)
      double precision dfault(*),tfault(*),diforc(*),wink(*),winkx(*),s(*)
      double precision stemp(*)
      character ofile*(*),pfile*(*)
c
c...  included dimension statements
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
c...  local variables
c
      integer ii,igroup,naxstp,nfirst
cdebug      integer idb
      double precision time
      logical*4 fulout,skc
c
c...  included variable definition statements
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
cdebug      write(6,*) "Hello from elastc_f!"
cdebug      write(6,*) "From elastc_f, ngauss, nsd, gauss:", ngauss,nsd,
cdebug     & (gauss(idb),idb=1,(nsd+1)*ngauss)
c
      if(idout.ne.0) open(kw,file=ofile,status="old",access="append")
      if(idsk.eq.0) open(kp,file=pfile,status="old",access="append")
      if(idsk.eq.1) open(kp,file=pfile,status="old",form="unformatted",
     & access="append")
      skc=iskopt.ge.0.and.iskopt.ne.1.and.numslp.ne.0
      call fill(bres,zero,neq)
      call fill(btot,zero,neq)
      call fill(gvec1,zero,neq)
      call fill(pvec,zero,neq)
      call fill(deld,zero,ndof*numnp)
      call fill(deldx,zero,ndof*numnp)
      call fill(dcur,zero,ndof*numnp)
      call fill(d,zero,ndof*numnp)
      call fill(dx,zero,ndof*numnp)
      call fill(dxcur,zero,ndof*numnp)
      call fill(deps,zero,nstr*ngauss*numel)
      if(numfn.ne.0) call fill(tfault,zero,numfn*ndof)
      if(ivisc.eq.1) call fill(beta,zero,nstr*ngauss*numel)
      if(iplas.eq.1) call fill(betb,zero,nstr*ngauss*numel)
      if(ipstrs.eq.1) call fill(st0,zero,nstr*nppts*numel)
      if(skc) call skclear(idslp,skew,numsn,nskdim,numnp)
      if(skc) call skcomp(x,d,skew,idslp,ipslp,nsd,ndof,nskdim,npdim,
     & ipstrs,numsn,numnp,nstep,lgdefp,kto)
c
      write(kto,600)
c*      call flush(kto)
      fulout=.true.
      ireform=1
      igroup=1
      nstep=0
      ntimdat(1)=nstep
      naxstp=0
      nittot=0
      ntimdat(7)=nittot
      nrftot=0
      ntimdat(8)=nrftot
      ndtot=0
      ntimdat(9)=ndtot
      ntimdat(10)=ireform
      call const(maxstp,delt,alfa,maxit,maxitc,lgdef,ibbar,utol,ftol,
     & etol,itmax,nintg,igroup,naxstp,nfirst,rtimdat,deltp,alfap,
     & ntimdat,nstep,maxitp,maxitcp,lgdefp,ibbarp,itmaxp,gtol)
c
c...transfer boundary conditions into global load vector btot(neq)
c   and displacement increment vector deld(ndof,numnp)
c
      call load(id,ibond,bond,dcur,deld,btot,histry,deltp,ndof,numnp,
     & neq,nhist,nstep,lastep,idout,kto,kw)
c
c...compute current split node displacements
c
      if(numfn.ne.0) then
        call loadf(fault,dfault,histry,deltp,nfault,nstep,ndof,numfn,
     &   nhist,lastep,idout,kto,kw)
      end if
c
c...add differential forces across internal free interfaces
c
      if(numdif.ne.0) call loadx(btot,diforc,histry,idx,idhist,ndof,neq,
     & numnp,nhist,nstep,lastep,idout,kto,kw)
c
c...compute forces due to prestresses, applied displacements, and split
c   nodes
c
      ii=1
      call elasmat(????????????)
      call formmat(stn,dmat,deps,beta,betb,prop,mat,iddmat,histry,
     & rtimdat,ntimdat,rgiter,ii,ngauss,nddmat,ndmat,
     & nprop,numat,ndof,nstr,numel,ipstrs,nhist,lastep,idout,kto,kw)
c************  still need to figure out how to do prestresses -- maybe
c************  with a separate program section.  Right now, it appears
c************  the best way is to compute the equivalent nodal forces and
c************  store these in a vector that gets subtracted each time.
c************  For now, ignore them and assume that they will be taken
c************  care of properly.
      if(nprestr.ne.0.and.ipstrs.eq.0) then
        call stresn(x,b,d,dx,tfault,stn,deps,beta,betb,scur,st0,dbeta,
     &   dbetb,skew,ien,lm,lmx,lmf,dmat,mat,prop,histry,infin,gauss,
     &   rtimdat,stol,iddmat,nen,numel,ndof,nsd,numnp,neq,nee,
     &   nstr,ngauss,nppts,ngem,nskdim,nhist,nprop,numat,numfn,numslp,
     &   numrot,lastep,nstep,lgdefp,ibbarp,ivisc,iplas,imhist,ipstrs,
     &   nprestr,nddmat,ndmat,idebug,idout,kto,kw,fulout)
      end if
c*************  Note that all routines need to be modified.  This shouldn't
c*************  be too hard in most cases.
      call formdf(b,x,d,dx,tfault,dcur,skew,prop,dmat,scur,histry,s,
     & stemp,gauss,ien,infin,lm,lmx,lmf,mat,iddmat,ngauss,nddmat,ndmat,
     & nprop,numat,nsd,ndof,nstr,nen,nee,neq,numel,numnp,numfn,numslp,
     & numrot,nskdim,ipstrs,nstep,lgdefp,ibbarp,nhist,lastep,idout,kto,
     & kw,ivisc,iplas,imhist)
      if(numfn.ne.0) call formf(b,x,d,dx,skew,histry,
     & ien,lm,lmx,lmf,gauss,
     & mat,infin,prop,dmat,stn,
     & nfault,dfault,tfault,
     & s,stemp,iddmat,
     & ngauss,nddmat,ndmat,nprop,numat,nsd,ndof,nstr,nen,nee,neq,numel,
     & numnp,numfn,numslp,numrot,nskdim,ipstrs,nstep,lgdefp,ibbarp,
     & nhist,lastep,idout,kto,kw,ivisc,iplas,imhist)
c
c...compute iterative solution to the elastic problem
c   return if convergence achieved or maximum iterations exceeded
c
c*************  The call to iterate needs to be updated for new element
c*************  methods.
c*************  Also note:  iterate could probably be called with
c*************  appropriate subroutine names to specify:
c     elastic solution
c     small/large strain
c     bbar or not
c     ?????
      call iterate(alnz,pcg,zcg,ja,                                     ! sparse
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
c...if ipstrs=1, equate prestresses with current stresses, and set
c   displacements to zero
c
      if(nprestr.ne.0.and.ipstrs.eq.1) then
c*        call fill(d,zero,ndof*numnp)
c*        call fill(dx,zero,ndof*numnp)
c*        call fill(deld,zero,ndof*numnp)
c*        call fill(deldx,zero,ndof*numnp)
c*        call fill(dcur,zero,ndof*numnp)
c*        call fill(dxcur,zero,ndof*numnp)
c*        call fill(deps,zero,nstr*ngauss*numel)
        call eqstrsql(st0,stn,x,gauss,ien,infin,nstr,ngauss,nppts,numel,
     &   nen,nsd,numnp,idout,kto,kw)
      end if
c
c...print all displacements, including faulted and slippery nodes
c
      time=zero
      if(idsk.eq.0) write(kp,700) nstep
      if(idsk.eq.0) write(kp,'(e15.4)') time
      if(idsk.eq.1) write(kp) nstep
      if(idsk.eq.1) write(kp) time
      call printd(d,deld,deltp,idslp,ndof,numnp,numnp,ione,idout,idsk,
     & kto,kw,kp)
      call printf(tfault,dfault,deltp,nfault,ndof,numfn,idout,idsk,kw,
     & kp)
      call printd(dx,deldx,deltp,idslp,ndof,numnp,numsn,itwo,idout,idsk,
     & kto,kw,kp)
c
c...print array telling whether each slippery node is locked
c   or free for the current time step
c
      call printl(idx,iwinkx,idslp,histry,ndof,numsn,numnp,nstep,nhist,
     & nwinkx,lastep,idsk,kp)
c
c...print the stresses and strains
c
      call prints(stn,eps,deps,beta,dbeta,betb,dbetb,nstr,ngauss,numel,
     & nstep,idout,idsk,kw,kp,ivisc,iplas)
      if(nintg.eq.1) then
        write(kto,800) ntimdat(7),ntimdat(8),ntimdat(9)
        if(idout.gt.0) write(kw,800) ntimdat(7),ntimdat(8),ntimdat(9)
      end if
      if(idout.ne.0) close(kw)
      close(kp)
c
600   format(//,'Beginning elastic solution:',/)
700   format('STEP ',i5)
800   format(/," Total number of equilibrium iterations        = ",i7,/,
     &         " Total number of stiffness matrix reformations = ",i7,/,
     &         " Total number of displacement subiterations    = ",i7)
      return
      end
c
c version
c $Id: elastc-tmp.f,v 1.1 2004/07/07 14:22:00 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
