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
      subroutine formes(x,d,dx,tfault,dmat,stn,skew,s,stemp,prop,gauss,
     & ien,lmx,lmf,iddmat,infin,n,ngauss,nddmat,nprop,nsd,ndof,nstr,nen,
     & nee,numnp,numfn,numslp,numrot,nskdim,iopt,ibbar,lgdef,idout,kto,
     & kw)
c
c...  subroutine to form the elemental stiffness matrix
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer infin,n,ngauss,nddmat,nprop,nsd,ndof,nstr,nen,nee,numnp
      integer numfn,numslp,numrot,nskdim,iopt,ibbar,lgdef,idout,kto,kw
      integer ien(nen),lmx(ndof,nen),lmf(nen)
      double precision x(nsd,numnp),d(ndof,numnp),dx(ndof,numnp)
      double precision tfault(ndof,numfn),dmat(nddmat,ngauss)
      double precision stn(nstr,ngauss),skew(nskdim,numnp),s(nee*nee)
      double precision stemp(nee*nee),prop(nprop),gauss(nsd+1,ngauss)
c
c...  included dimension and typing statements
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
      double precision xl(24),dl(24)
c
cdebug      write(6,*) "Hello from formes_f!"
cdebug      write(6,*) "From formes_f, ngauss,nsd,gauss: ",ngauss,nsd,
cdebug     & ((gauss(idb,jdb),idb=1,nsd+1),jdb=1,ngauss)
c
      call fill(s,zero,nee*nee)
c
c...  localize coordinates and update them for large deformations, if
c     required
c
      call lcoord(x,xl,ien,nen,nsd,numnp)
      if(lgdef.gt.izero) call ldupdat(d,dx,tfault,dl,xl,ien,lmx,lmf,
     & ndof,nsd,nen,numnp,numfn,numslp,iopt,lgdef)
c
c...  construct local stiffness matrix, symmetrize it, and rotate for
c     skew boundary conditions
c
      call stiffql(s,stn,dmat,xl,prop,gauss,ien,iddmat,infin,n,ngauss,
     & nddmat,nprop,ndof,nsd,nstr,nen,nee,ibbar,lgdef,idout,kto,kw)
      if(numrot.ne.izero) call rstiff(s,stemp,skew,ien,ndof,numnp,nen,
     & nee,nskdim)
      return
      end
c
c version
c $Id: formes.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
