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
      subroutine esfcomp(stn,stnc,ee,beta,dbeta,betb,dbetb,prop,
     & rtimdat,stol,n,nstr,ndof,nprop,ipstrs,nstep,lgdefp,idout,kto,kw,
     & ivisc,iplas)
c
c...subroutine to compute the stresses corresponding to the zero of
c   the effective stress function (ref. Kojic and Bathe, 1987, Int.
c   J. Num. Meth. Eng., v.24, pp 1509-1532).
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer n,nstr,ndof,nprop,ipstrs,nstep,lgdefp,idout,kto,kw
      integer ivisc,iplas
      double precision stn(nstr),stnc(nstr),ee(nstr),beta(nstr)
      double precision dbeta(nstr),betb(nstr),dbetb(nstr),prop(nprop)
      double precision stol
c
c...  included dimension and type statements
c
      include "rtimdat_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c... intrinsic functions
c
      intrinsic log,sqrt
c
c...  user-defined functions
c
      double precision sprod,rtsafe
      external sprod,rtsafe
c
c...  local variables
c
      integer i,iopd
      double precision brac,e,pois,emhu,anpwr,alpha,rk0,hards,r1,r2,r3
      double precision r3a,r4,r5
      double precision ae,be,yield,eprstrn,efsts,emeant,emeanp,emeantp
      double precision emeantpp,efstsi,ds,d,efstse,comp,dlam,gam,strtau
      double precision rm,dl1,dl2,t1,t2,bs,c,add,flo,fhi,smean,dnm,g1,g2
      double precision g3,g4,sdev,stau,sinv1
      double precision devsts(6),eep(6),betat(6),betbt(6)
      logical succes,test,visc,plas
c
c...  included variable definitions
c
      include "rtimdat_def.inc"
c
cdebug2      write(6,*) "Hello from esfcomp_f!"
cdebug2      write(6,*) n,nstr,ndof,nprop,ipstrs,nstep,lgdefp,idout,kto,kw,
cdebug2     & ivisc,iplas
      test=.false.
      visc=.false.
      plas=.false.
      if(ivisc.eq.1) visc=.true.
      if(iplas.eq.1) plas=.true.
      brac=5000.0d0
      e=prop(1)
      pois=prop(2)
      if(ipstrs.eq.1.and.nstep.eq.izero.and.tpois.gt.zero) pois=tpois
      emhu=prop(4)
      anpwr=prop(5)
      alpha=prop(6)
      rk0=prop(7)
      hards=prop(8)
      r1=alpha+one/root3
      r2=one-two*pois
      r3=three*alpha*e
      r3a=r3*alpha
      r4=three*alpha*alpha+half
      r5=one+pois
      ae=r5/e
      be=r2/e
      yield=rk0/r1
      call fill(betat,zero,nstr)
      call fill(betbt,zero,nstr)
      if(visc.and.nstep.gt.izero) call dcopy(nstr,beta,ione,betat,ione)
      if(plas.and.nstep.gt.izero) call dcopy(nstr,betb,ione,betbt,ione)
      if(lgdefp.gt.izero) then
        do i=1,ndof
          betat(i)=log(one+betat(i))
          betbt(i)=log(one+betbt(i))
        end do
      end if
      if(plas.and.hards.ge.zero) yield=yield+hards*r1*
     & sqrt(two*sprod(betbt,betbt)/r4)
cdebug2      write(6,*) "Point 1"
c
c...compute deviatoric strain and (initial) effective stress
c
      call invar(devsts,sinv1,efsts,stn)
      eprstrn=prop(nprop-1)
      if(lgdefp.gt.izero.and.eprstrn.gt.-one) eprstrn=log(one+eprstrn)
      emeant=third*(ee(1)+ee(2)+ee(3))
      emeanp=zero
      if(plas) emeanp=third*(betbt(1)+betbt(2)+betbt(3))
      emeantp=emeant-emeanp
      emeantpp=emeantp-third*eprstrn
      efstsi=efsts
cdebug2      write(6,*) "Point 2"
      do i=1,nstr
        eep(i)=ee(i)-betat(i)-betbt(i)
        if(i.le.3) eep(i)=eep(i)-emeantp
      end do
      ds=sprod(eep,eep)
cdebug2      write(6,*) "Point 3"
      d=sqrt(ds)
      efstse=d/ae
c
c...  compute effective stress, plastic strain, and mean stress for
c     different cases.  Use effective stress function, if necessary
c
      if(hards.le.zero.and.(.not.visc).and.alpha.eq.zero) then
        comp=root3*efstse
        test=.false.
        if(comp.ge.yield.and.plas) then
          efsts=yield/root3
          dlam=two*(d-ae*efsts)/root3
          gam=zero
        else
          efsts=efstse
          dlam=zero
          gam=zero
        end if
cdebug2        write(6,*) "Point 4"
      else if(hards.le.zero.and.visc.and.alpha.eq.zero) then
        comp=root3*efstse
        if(comp.ge.yield.and.plas) then
          efsts=yield/root3
          strtau=(one-alfap)*efstsi+alfap*efsts
          gam=half*(strtau/emhu)**(anpwr-one)/emhu
          rm=deltp*gam*(one-alfap)
          dlam=two*(sqrt(ds-two*rm*sprod(eep,devsts)-rm*rm*efstsi*
     &     efstsi)-efsts*(ae-deltp*gam*alfap))/root3
          test=.false.
        else
          test=.true.
          efsts=efstse
          dl1=zero
          dl2=zero
        end if
cdebug2        write(6,*) "Point 5"
      else if(hards.gt.zero.and.(.not.visc).and.alpha.eq.zero) then
        comp=root3*efstse
        test=.false.
        if(comp.ge.yield.and.plas) then
          efsts=(root3*yield+two*hards*d)/(two*ae*hards+three)
          dlam=(root3*efsts-yield)/hards
          gam=zero
        else
          efsts=efstse
          dlam=zero
          gam=zero
        end if
cdebug2        write(6,*) "Point 6"
      else if(hards.gt.zero.and.visc.and.alpha.eq.zero) then
        test=.true.
        efsts=efstse
        dl1=root3/hards
        dl2=-(yield/hards)
cdebug2        write(6,*) "Point 7"
      else if(hards.le.zero.and.(.not.visc).and.alpha.gt.zero) then
        dlam=r1*emeantpp/alpha+r1*r2*(efstse-r1*yield)/r3a
        test=.false.
        if(dlam.ge.zero.and.plas) then
          efsts=(r3*(two*alpha*d+emeantpp)+r1*r2*yield)/
     &     (six*alpha*alpha*r5+r2)
          dlam=two*r1*(d-ae*efsts)
          gam=zero
        else
          efsts=efstse
          dlam=zero
          gam=zero
        end if
cdebug2        write(6,*) "Point 8"
      else if(hards.le.zero.and.visc.and.alpha.gt.zero) then
        test=.true.
        efsts=efstse
        dl1=r1*r2/r3a
        dl2=r1*emeantpp/alpha-dl1*r1*yield
cdebug2        write(6,*) "Point 9"
      else if(hards.gt.zero.and.(.not.visc).and.alpha.gt.zero) then
        test=.false.
        dlam=r1*(r3*emeantpp+r2*(efstse-r1*yield))/(hards*r2*r1*r1+r3a)
        if(dlam.ge.zero.and.plas) then
          rm=hards*r2*r1*r1+r3a
          efsts=(two*d*rm-r3*emeantpp-r2*r1*yield)/(two*ae*rm+r2)
          dlam=two*r1*(d-ae*efsts)
          gam=zero
        else
          efsts=efstse
          dlam=zero
          gam=zero
        end if
cdebug2        write(6,*) "Point 10"
      else if(hards.gt.zero.and.visc.and.alpha.gt.zero) then
        test=.true.
        efsts=efstse
        rm=hards*r2*r1*r1+r3a
        dl1=r1*r2/rm
        dl2=r1*(r3*emeantpp-r2*r1*yield)/rm
cdebug2        write(6,*) "Point 11"
      end if
c
c...  use effective stress function, if required
c
      if(test) then
        t1=two*r1
        t2=deltp*alfap
        bs=deltp*deltp*(one-alfap)*(one-alfap)*efstsi*efstsi
        c=two*deltp*(one-alfap)*sprod(eep,devsts)
        if(efsts.lt.efstsi) efsts=efstsi
        iopd=1
        add=brac*stol*efsts
        if(add.eq.zero) add=stol
        flo=zero
        fhi=efsts+add
        if(.not.plas) then
          dl1=zero
          dl2=zero
        end if
cdebug2        write(6,*) "Point 12"
c
c...bracket the root
c
        call zbrac(flo,fhi,succes,ae,bs,c,ds,dl1,dl2,t1,t2,emhu,
     &   anpwr,efstsi,gam,dlam,deltp,alfap,iopd,plas)
cdebug2        write(6,*) "Point 13"
c
c...compute effective stress
c
        if(succes) then
          iopd=2
          efsts=rtsafe(flo,fhi,stol,ae,bs,c,ds,dl1,dl2,t1,t2,emhu,
     &     anpwr,efstsi,gam,dlam,deltp,alfap,iopd,idout,kto,kw,plas)
cdebug2          write(6,*) "Point 14"
        else
          write(kto,800) n,flo,fhi,ae,bs,c,ds,dl1,dl2,t1,t2,emhu,anpwr,
     &     efstsi,gam,dlam
cdebug2          write(6,*) "Point 15"
          stop
        end if
      end if
c
c...Compute increment of creep strain, plastic strain, deviatoric
c   stress, and stress
c
      emeantpp=emeantpp-alpha*dlam/r1
      smean=emeantpp/be
      dnm=ae+alfap*deltp*gam
      if(efsts.ne.zero) dnm=dnm+dlam/(two*r1*efsts)
      g1=one/dnm
      g2=(one-alfap)*deltp*gam
      g3=deltp*gam
      g4=dlam/r1
      do i=1,nstr
        sdev=g1*(eep(i)-g2*devsts(i))
        stnc(i)=sdev
        if(i.le.3) stnc(i)=smean+sdev
        stau=(one-alfap)*devsts(i)+alfap*sdev
        if(visc) dbeta(i)=g3*stau
        if(plas) dbetb(i)=zero
        if(plas.and.efsts.ne.zero) dbetb(i)=g4*sdev/(two*efsts)
        if(plas.and.i.le.3) dbetb(i)=dbetb(i)+g4*alpha
      end do
cdebug2      write(6,*) "Point 16"
800   format('Root not bracketed for element # ',i7,/,
     & ' flo   = ',1pe15.8,' fhi   = ',1pe15.8,' ae    = ',1pe15.8,/,
     & ' bs    = ',1pe15.8,' c     = ',1pe15.8,' ds    = ',1pe15.8,/,
     & ' dl1   = ',1pe15.8,' dl2   = ',1pe15.8,' t1    = ',1pe15.8,/,
     & ' t2    = ',1pe15.8,' emhu  = ',1pe15.8,' anpwr = ',1pe15.8,/,
     & ' efstsi= ',1pe15.8,' gam   = ',1pe15.8,' dlam  = ',1pe15.8)
      return
      end
c
c version
c $Id: esfcomp.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
