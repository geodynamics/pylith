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
      subroutine elas_strs_cmp_ss(
     & b,neq,                                                           ! force
     & x,d,dx,nsd,ndof,numnp,                                           ! global
     & tfault,numfn,numslp,                                             ! fault
     & state,dmat,ien,lm,lmx,lmf,infiel,nstr,nddmat,nstatesz,ndmatsz,   ! elemnt
     & numelt,nconsz,                                                   ! elemnt
     & infmat,matgpt,elas_strs,                                         ! materl
     & gauss,sh,shj,infetype,netypes,                                   ! eltype
     & skew,nskdim,numrot,                                              ! skew
     & getshape,bmatrix,idebug,idout,kto,kw,ierr)                       ! info
c
c...program to compute the total stress and strain for the current
c   iteration for a given material model
c
      include "implicit.inc"
c
c...  dimension parameters
c
      include "nshape.inc"
c
c...  subroutine arguments
c
      integer neq,nsd,ndof,numnp,numfn,numslp,nstr,nddmat,nstatesz
      integer ndmatsz,numelt,nconsz,matgpt,netypes,nskdim,numrot,idebug
      integer idout,kto,kw,ierr
      integer ien(nconsz),lm(ndof,nconsz),lmx(ndof,nconsz),lmf(nconsz)
      integer infiel(6,numelt),infmat(6),infetype(4,netypes)
      double precision b(neq),x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn)
      double precision state(nstr,nstatesz),dmat(nddmat,ndmatsz)
      double precision gauss(nsd+1,ngaussmax,netypes)
      double precision sh(nsd+1,nenmax,ngaussmax,netypes)
      double precision shj(nsd+1,nenmax,ngaussmax,netypes)
      double precision skew(nskdim,numnp)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  external routines
c
      external elas_strs,getshape,bmatrix
c
c...  local variables
c
      integer npage,nmatel,nstate,ind,iel,indien,ietype,indstate,inddmat
      integer ngauss,nen,nee,l,indstateg,inddmatg
      double precision dl(60),xl(60),scur(162),ee(162),p(60),det(27)
      logical debug
c
cdebug      write(6,*) "Hello from elas_strs_cmp_ss_f!"
c
      debug=(idebug.eq.1).and.(idout.gt.1)
      npage=50
      nmatel=infmat(2)
      nstate=infmat(3)
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
c
c...  localize coordinates and displacements
c
        call lcoord(x,xl,ien(indien),nen,nsd,numnp)
        call ldisp(dl,d,ien(indien),ndof,nen,numnp)
        if(numfn.ne.0) call adfldp(dl,lmf(indien),tfault,ndof,nen,numfn)
        if(numslp.ne.0) call addsn(dl,dx,ien,lmx(1,indien),ndof,nen,
     &   numnp)
c
c...  compute strains
c
        call bdeld_ss(xl,dl,sh(1,1,1,ietype),shj(1,1,1,ietype),ee,det,
     &   gauss(1,1,ietype),iel,nen,nee,nsd,ndof,nstr,ngauss,ierr,
     &   getshape,bmatrix)
        if(ierr.ne.0) return
c
c...  loop over gauss points, compute stresses, and transfer them into
c     scur
c
        do l=1,ngauss
          indstateg=indstate+(l-1)*nstate
          inddmatg=inddmat+(l-1)*nddmat
          call elas_strs(state(1,indstateg),ee(nstr*(l-1)),
     &     dmat(1,inddmatg),nstr,nstate,nddmat)
          call dcopy(nstr,state(1,indstateg),ione,scur(nstr*(l-1)),ione)
        end do
c
c...  compute equivalent nodal loads
c
        call fill(p,zero,ndof*nen)
        call eforce(xl,sh(1,1,1,ietype),shj(1,1,1,ietype),det,
     &   gauss(1,1,ietype),scur,p,iel,nen,nsd,ndof,nstr,ngauss,ierr,
     &   getshape,bmatrix)
        if(ierr.ne.0) return
        if(numrot.ne.0) call rpforc(p,skew,ien(indien),ndof,numnp,nen,
     &   nskdim)
        if(debug) then
          if(ind.eq.1.or.mod(ind,npage).eq.0) write(kw,1000)
          call prntforc(iel,p,ien(indien),nen,ndof,idout,kw)
        end if
        call addfor(b,p,lm(1,indien),lmx(1,indien),neq,nee)
      end do
 1000 format(//," local forces computed from stress field",//)
      return
      end
c
c version
c $Id: elas_strs_cmp_ss.f,v 1.1 2004/06/16 16:53:36 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
