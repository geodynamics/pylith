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
      subroutine elas_strs_mat_cmp_ss(
     & alnz,ja,nnz,                                                     ! sparse
     & bintern,neq,                                                     ! force
     & x,d,numnp,                                                       ! global
     & dx,numslp,numsn,                                                 ! slip
     & tfault,numfn,                                                    ! fault
     & s,stemp,                                                         ! stiff
     & state,dstate,state0,dmat,ien,lm,lmx,lmf,infiel,iddmat,nstatesz,  ! elemnt
     & nstatesz0,ndmatsz,numelt,nconsz,nprestrflag,ipstrs,ipauto,       ! elemnt
     & nstate0,                                                         ! elemnt
     & prop,nmatel,imatvar,nstate,nprop,matgpt,elas_strs_mat,           ! materl
     & td_strs_mat,matchg,tminmax,                                      ! materl
     & gauss,sh,shj,infetype,                                           ! eltype
     & rtimdat,ntimdat,rgiter,                                          ! timdat
     & skew,numrot,                                                     ! skew
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...  subroutine to compute the total stress and strain for the current
c     iteration for a given material model and update the global load
c     vector and the material and stiffness matrices.
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
      integer nnz,neq,numnp,numslp,numsn,numfn,nstatesz,nstatesz0
      integer ndmatsz,numelt,nconsz,nprestrflag,ipstrs,ipauto,nstate0
      integer nmatel,imatvar,nstate,nprop,matgpt,numrot,ierr
      integer ja(nnz),ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz)
      integer lmf(nconsz),infiel(7,numelt),iddmat(nstr,nstr)
      integer infetype(4,netypes)
      character errstrng*(*)
      logical matchg
      double precision alnz(nnz),bintern(neq),x(nsd,numnp),d(ndof,numnp)
      double precision dx(ndof,numnp),tfault(ndof,numfn)
      double precision s(neemax*neemax),stemp(neemax*neemax)
      double precision state(nstr,nstatesz),dstate(nstr,nstatesz)
      double precision state0(nstate0,nstatesz0)
      double precision dmat(nddmat,ndmatsz),prop(nprop),tminmax
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
c
c...  external routines
c
      external elas_strs_mat,td_strs_mat,getshape,bmatrix
c
c...  local variables
c
      integer ind,iel,indien,ietype,indstate,inddmat,indstate0
      integer ngauss,nen,nee,l,indstateg,inddmatg,indstate0g,ngtest
      integer ngaussdim
      double precision dl(60),xl(60),scur(162),ee(162),p(60),det(27)
c
cdebug      integer idb,jdb
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
      include "rgiter_def.inc"
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "Hello from elas_strs_mat_cmp_ss_f!"
c
      ngtest=izero
      if(imatvar.eq.izero) ngtest=ngaussmax
c
c...  loop over elements in a material group
c
      do ind=matgpt,matgpt+nmatel-1
        iel=infiel(4,ind)
        indien=infiel(1,iel)
        ietype=infiel(3,iel)
        indstate=infiel(5,iel)
        inddmat=infiel(6,iel)
        indstate0=infiel(7,iel)
        ngauss=infetype(1,ietype)
        nen=infetype(2,ietype)
        nee=infetype(4,ietype)
        indstateg=indstate
        inddmatg=inddmat
        indstate0g=indstate0
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
cdebug          write(6,*) "ind,iel,l,indstate0g,state0:",ind,iel,l,
cdebug     &     indstate0g,(state0(idb,indstate0g),idb=1,nstate0)
          call elas_strs_mat(dstate(1,indstateg),state0(1,indstate0g),
     &     ee(nstr*(l-1)+1),dmat(1,inddmatg),nstate,nstate0,ierr,
     &     errstrng)
          if(ierr.ne.izero) return
          call dcopy(nstr,dstate(1,indstateg),ione,scur(nstr*(l-1)+1),
     &     ione)
          indstateg=indstateg+nstate
          inddmatg=inddmatg+ione
          indstate0g=indstate0g+ione
        end do
c
c...  compute equivalent nodal loads and add them to global load
c
        call fill(p,zero,ndof*nen)
        call eforce(xl,sh(1,1,1,ietype),shj(1,1,1,ietype),det,
     &   gauss(1,1,ietype),scur,p,iel,nen,ngauss,getshape,bmatrix,ierr,
     &   errstrng)
        if(ierr.ne.izero) return
        if(numrot.ne.izero) call rpforc(p,skew,ien(indien),numnp,nen)
        call addfor(bintern,p,lm(1,indien),lmx(1,indien),neq,nee)
c
c...  compute element stiffness and add it to global stiffness
c
        ngaussdim=max(ngauss,ngtest)
        call formes_ss(
     &   x,numnp,                                                       ! global
     &   s,stemp,                                                       ! stiff
     &   dmat(1,inddmat),ien(indien),lm(1,indien),iddmat,iel,           ! elemnt
     &   gauss(1,1,ietype),sh(1,1,1,ietype),shj(1,1,1,ietype),          ! eltype
     &   ngauss,ngaussdim,nen,nee,                                      ! eltype
     &   skew,numrot,                                                   ! skew
     &   getshape,bmatrix,                                              ! bbar
     &   ierr,errstrng)                                                 ! errcod
        if(ierr.ne.izero) return
        call addstf(alnz,s,lm(1,indien),lmx(1,indien),ja,nee,numsn,nnz)
      end do
      return
      end
c
c version
c $Id: elas_strs_mat_cmp_ss.f,v 1.10 2005/02/23 23:54:56 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
