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
      subroutine td_strs_cmp_ss(
     & b,neq,                                                           ! force
     & x,d,numnp,                                                       ! global
     & dx,numslp,                                                       ! slip
     & tfault,numfn,                                                    ! fault
     & state,dstate,dmat,ien,lm,lmx,lmf,infiel,nstatesz,ndmatsz,        ! elemnt
     & numelt,nconsz,                                                   ! elemnt
     & prop,nmatel,nstate,nprop,matgpt,elas_strs,td_strs,               ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & rtimdat,ntimdat,rgiter,                                          ! timdat
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...program to compute the total stress and strain for the current
c   iteration for a given material model
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer neq,numnp,numslp,numfn,nstatesz,ndmatsz,numelt,nconsz
      integer nmatel,nstate,nprop,matgpt,numrot,ierr
      integer ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz),lmf(nconsz)
      integer infiel(6,numelt),infetype(4,netypes)
      character errstrng*(*)
      double precision b(neq),x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn)
      double precision state(nstr,nstatesz),dstate(nstr,nstatesz)
      double precision dmat(nddmat,ndmatsz),prop(nprop)
      double precision gauss(nsd+1,ngaussmax,netypes)
      double precision sh(nsd+1,nenmax,ngaussmax,netypes)
      double precision shj(nsd+1,nenmax,ngaussmax,netypes)
      double precision skew(nskdim,numnp)
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
      include "ntimdat_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  external routines
c
      external elas_strs,td_strs,getshape,bmatrix
c
c...  local variables
c
      integer ind,iel,indien,ietype,indstate,inddmat
      integer ngauss,nen,nee,l,indstateg,inddmatg,incstate
      double precision dl(60),xl(60),scur(162),ee(162),p(60),det(27)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from td_strs_cmp_ss_f!"
c
      incstate=nstr*nstate
c
c...  loop over elements in a material group
c
      do ind=matgpt,matgpt+nmatel-1
        iel=infiel(4,ind)
        indien=infiel(1,iel)
        ietype=infiel(3,iel)
        indstate=infiel(5,iel)
        inddmat=infiel(6,iel)
        ngauss=infetype(1,ietype)
        nen=infetype(2,ietype)
        nee=infetype(4,ietype)
        indstateg=indstate
        inddmatg=inddmat
c
c...  localize coordinates and displacements
c
        call lcoord(x,xl,ien(indien),nen,numnp)
        call ldisp(dl,d,ien(indien),nen,numnp)
        if(numfn.ne.0) call adfldp(dl,lmf(indien),tfault,nen,numfn)
        if(numslp.ne.0) call addsn(dl,dx,ien,lmx(1,indien),nen,numnp)
c
c...  compute strains
c
        call bdeld_ss(xl,dl,sh(1,1,1,ietype),shj(1,1,1,ietype),ee,det,
     &   gauss(1,1,ietype),iel,nen,nee,ngauss,getshape,bmatrix,ierr,
     &   errstrng)
        if(ierr.ne.izero) return
c
c...  loop over gauss points, compute stresses, and transfer them into
c     scur
c
        do l=1,ngauss
          call td_strs(state(1,indstateg),dstate(1,indstateg),
     &     ee(nstr*(l-1)),dmat(1,inddmatg),prop,rtimdat,rgiter,ntimdat,
     &     nstate)
          call dcopy(nstr,dstate(1,indstateg),ione,scur(nstr*(l-1)),
     &     ione)
          indstateg=indstateg+incstate
          inddmatg=inddmatg+nddmat
        end do
c
c...  compute equivalent nodal loads
c
        call fill(p,zero,ndof*nen)
        call eforce(xl,sh(1,1,1,ietype),shj(1,1,1,ietype),det,
     &   gauss(1,1,ietype),scur,p,iel,nen,ngauss,getshape,bmatrix,ierr,
     &   errstrng)
        if(ierr.ne.izero) return
        if(numrot.ne.izero) call rpforc(p,skew,ien(indien),numnp,nen)
        call addfor(b,p,lm(1,indien),lmx(1,indien),neq,nee)
      end do
      return
      end
c
c version
c $Id: td_strs_cmp_ss.f,v 1.3 2004/07/09 01:33:38 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
