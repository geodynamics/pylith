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
      subroutine formmat(
     & stn,dmat,eps,beta,betb,                                          ! stress
     & prop,mat,iddmat,                                                 ! elemnt
     & histry,rtimdat,ntimdat,rgiter,iter,                              ! timdat
     & ngauss,nddmat,ndmat,nprop,numat,ndof,nstr,numel,ipstrs,nhist,    ! dimens
     & lastep,idout,kto,kw)                                             ! dimens
c
c...subroutine to reform the d-matrix
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer iter,ngauss,nddmat,ndmat,nprop,numat,ndof,nstr,numel
      integer ipstrs,nhist,lastep,idout,kto,kw
      integer mat(numel)
      double precision stn(nstr,ngauss,numel),dmat(nddmat,ngauss,ndmat)
      double precision eps(nstr,ngauss,numel),beta(nstr,ngauss,numel)
      double precision betb(nstr,ngauss,numel),prop(nprop,numat)
      double precision histry(nhist,lastep+1)
c
c...  included dimension and type statements
c
      include "iddmat_dim.inc"
      include "ntimdat_dim.inc"
      include "rtimdat_dim.inc"
      include "rgiter_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
c
c...  local variables
c
cdebug      integer idb,jdb,kdb
      integer m,l,n
      double precision ptmp(30)
c
c...  included variable definitions
c
      include "ntimdat_def.inc"
c
cdebug      write(6,*) "From formmat_f, ngauss: ",ngauss
cdebug      write(6,*) "Hello from formmat_f!"
cdebug2      write(6,*) "iter,ngauss,nddmat,ndmat,nprop,numat,ndof,nstr,numel"
cdebug2      write(6,*) "ipstrs,nhist,lastep,idout,kto,kw"
cdebug2      write(6,*) iter,ngauss,nddmat,ndmat,nprop,numat,ndof,nstr,numel
cdebug2      write(6,*) ipstrs,nhist,lastep,idout,kto,kw
cdebug2      write(6,*) "stn:",(((stn(idb,jdb,kdb),idb=1,nstr),jdb=1,ngauss),
cdebug2     & kdb=1,numel)
cdebug2      write(6,*) "dmat:",(((dmat(idb,jdb,kdb),idb=1,nddmat),
cdebug2     & jdb=1,ngauss),kdb=1,ndmat)
cdebug2      write(6,*) "eps:",(((eps(idb,jdb,kdb),idb=1,nstr),jdb=1,ngauss),
cdebug2     & kdb=1,numel)
cdebug2      write(6,*) "beta:",(((beta(idb,jdb,kdb),idb=1,nstr),jdb=1,ngauss),
cdebug2     & kdb=1,numel)
cdebug2      write(6,*) "prop:",((prop(idb,jdb),idb=1,nprop),jdb=1,numat)
cdebug2      write(6,*) "histry:",((histry(idb,jdb),idb=1,nhist),
cdebug2     & jdb=1,lastep+1)
c
      if((ivisc.eq.izero).and.(iplas.eq.izero).and.iter.eq.ione) then
        do m=1,numat
          call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &     idout,kto,kw,imhist)
          do l=1,ngauss
            call matinit(stn(1,l,m),eps(1,l,m),beta(1,l,m),betb(1,l,m),
     &       dmat(1,l,m),ptmp,rtimdat,iddmat,nstr,nddmat,nprop,
     &       ndof,ipstrs,nstep,lgdefp,ivisc,iplas)
cdebug            write(6,"(2i7,21(2x,1pe15.8))") m,l,(dmat(jdb,l,m),jdb=1,nddmat)
          end do
        end do
      else if(((ivisc.eq.ione).or.(iplas.eq.ione)).and.
     & (iter.eq.ione.or.(iter.gt.ione.and.nstep.gt.izero))) then
        do n=1,numel
          m=mat(n)
          call mathist(ptmp,prop(1,m),histry,nprop,m,nstep,nhist,lastep,
     &     idout,kto,kw,imhist)
          do l=1,ngauss
            if(iter.eq.ione) call matinit(stn(1,l,n),eps(1,l,n),
     &       beta(1,l,n),betb(1,l,n),dmat(1,l,n),ptmp,rtimdat,
     &       iddmat,nstr,nddmat,nprop,ndof,ipstrs,nstep,lgdefp,ivisc,
     &       iplas)
            if(iter.gt.ione) call matprtb(stn(1,l,n),eps(1,l,n),
     &       beta(1,l,n),betb(1,l,n),dmat(1,l,n),ptmp,rtimdat,
     &       rgiter,iddmat,n,nstr,ndof,nprop,nddmat,ipstrs,nstep,lgdefp,
     &       idout,kto,kw,ivisc,iplas)
cdebug2            write(6,*) "n,l,dmat:",n,l,(dmat(jdb,l,n),jdb=1,nddmat)
          end do
        end do
      end if
      return
      end
c
c version
c $Id: formmat.f,v 1.1 2004/07/02 18:42:02 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
