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
      subroutine matinit(stn,eps,beta,betb,dmat,prop,rtimdat,
     & iddmat,nstr,nddmat,nprop,ndof,ipstrs,nstep,lgdefp,ivisc,iplas)
c
c...constructs the tangent material matrix relating stress,strain
c   this routine is used at the beginning of a time step
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,nddmat,nprop,ndof,ipstrs,nstep,lgdefp,ivisc,iplas
      double precision stn(nstr),eps(nstr),beta(nstr),betb(nstr)
      double precision dmat(nddmat),prop(nprop)
c
c...  included dimension and type statements
c
      include "iddmat_dim.inc"
      include "rtimdat_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  intrinsic functions
c
      intrinsic log,sqrt
c
c...  user-defined functions
c
      double precision sprod
      external sprod
c
c...  local variables
c
      integer i
      double precision e,pois,emhu,anpwr,alpha,rk0,hards,ae,b1,gam,dlam
      double precision sinv1,steff,r1,r2,etm,yield,eeffp,strs,emeanp
      double precision emeant,emeantp,dnm,cp,cm
      double precision etmp(6),sdev(6),betat(6),betbt(6)
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
c
c*      write(6,*) "Hello from matinit_f!"
c
      call fill(dmat,zero,nddmat)
      e=prop(1)
      pois=prop(2)
      if(ipstrs.eq.1.and.nstep.eq.0.and.tpois.gt.zero) pois=tpois
      emhu=prop(4)
      anpwr=prop(5)
      alpha=prop(6)
      rk0=prop(7)
      hards=prop(8)
      ae=(one+pois)/e
      b1=one-two*pois
      gam=zero
      dlam=zero
      if(nstep.gt.0.and.((ivisc.eq.1).or.(iplas.eq.1))) then
        call invar(sdev,sinv1,steff,stn)
        if(ivisc.eq.1) then
          gam=(steff/emhu)**(anpwr-one)
          gam=half*gam/emhu
        end if
        if(iplas.eq.1) then
          r1=alpha+one/root3
          r2=sqrt(three*alpha*alpha+half)
          etm=sqrt(two)*r1/r2
          yield=rk0/r1
          call fill(betat,zero,nstr)
          if(ivisc.eq.1) call dcopy(nstr,beta,ione,betat,ione)
          call dcopy(nstr,betb,ione,betbt,ione)
          if(lgdefp.gt.0) then
            do i=1,ndof
              betat(i)=log(one+betat(i))
              betbt(i)=log(one+betbt(i))
            end do
          end if
          eeffp=etm*sqrt(sprod(betbt,betbt))
          if(hards.gt.zero) yield=yield+hards*eeffp
          strs=(alpha*sinv1+steff)/r1
          if(strs.le.yield) then
            dlam=zero
          else 
            emeant=third*(eps(1)+eps(2)+eps(3))
            emeanp=third*(betbt(1)+betbt(2)+betbt(3))
            emeantp=emeant-emeanp
            do i=1,nstr
              etmp(i)=eps(i)-betat(i)-betbt(i)
              if(i.le.3) etmp(i)=etmp(i)-emeantp
            end do
            dlam=two*r1*sprod(etmp,etmp)-steff*r1*(ae+deltp*gam)
c******  DOUBLE CHECK THE ABOVE EQUATION!!!!!!
            if(dlam.lt.zero) dlam=zero
          end if
        end if
      end if
      dnm=ae+dlam+deltp*gam
      cp=one/dnm
      cm=e/b1
      dmat(iddmat(1,1))=third*(two*cp+cm)
      dmat(iddmat(1,2))=third*(cm-cp)
      dmat(iddmat(1,3))=dmat(iddmat(1,2))
      dmat(iddmat(2,2))=dmat(iddmat(1,1))
      dmat(iddmat(2,3))=dmat(iddmat(1,2))
      dmat(iddmat(3,3))=dmat(iddmat(1,1))
      dmat(iddmat(4,4))=half*cp
      dmat(iddmat(5,5))=dmat(iddmat(4,4))
      dmat(iddmat(6,6))=dmat(iddmat(4,4))
      return
      end
c
c version
c $Id: matinit.f,v 1.1 2004/07/07 15:44:31 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
