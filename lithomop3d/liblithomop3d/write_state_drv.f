c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
c  Copyright 2005 Rensselaer Polytechnic Institute.
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
      subroutine write_state_drv(
     & state,dstate,ivfamily,nvfamilies,numelv,nstatesz,                ! elemnt
     & infmatmod,                                                       ! materl
     & ngauss,                                                          ! eltype
     & delt,nstep,                                                      ! timdat
     & istatout,                                                        ! ioopts
     & idout,idsk,iucd,kw,kp,kucd,ucdroot,iprestress)                   ! ioinfo
c
c...  program to print state variables
c     Note:  at present, it is assumed that the same istatout array
c     is used for the elastic and time-dependent solutions.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "materials.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nvfamilies,numelv,nstatesz,nstep,idout,idsk,iucd,kw,kp
      integer ivfamily(5,nvfamilies),infmatmod(6,nmatmodmax)
      integer istatout(3,nstatesmax)
      double precision delt
      double precision state(nstatesz),dstate(nstatesz)
c
c...  included dimension and type statements
c
      include "labels_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer itmp1(nstatesmax,3),itmp2(3*nstatesmax)
      integer nstatestot,nstatesinc,nstatesrate,nstatesout
      integer ifam,nelfamily,matmodel,indstate,nstate
      integer matmodpt(nstatesmax,nmatmodmax)
      integer m50,i,iel,imat,ietype,ngauss,matmodel,nstate,indstate,l
      integer indstateg,nout,j
      double precision stmp(6),tmult
      character statedescr*21
cdebug      integer idb,jdb
c
c...  included variable definitions
c
      include "labels_def.inc"
c
cdebug      write(6,*) "Hello from write_state_f!"
c
      if(delt.gt.zero) then
        delti=one/delt
      else
        delti=one
      end if
c
c...  create index array of state variables to output
c
      call ifill(itmp1,izero,3*nstatesmax)
      call ifill(itmp2,izero,3*nstatesmax)
      nstatestot=izero
      nstatesinc=izero
      nstatesrate=izero
      nstatesout=izero
      do i=1,nstatesmax
        if(istatout(1,i).ne.izero) then
          nstatestot=nstatestot+ione
          itmp(nstatestot,1)=i
        end if
        if(istatout(2,i).ne.izero) then
          nstatesinc=nstatesinc+ione
          itmp(nstatesinc,2)=i
        end if
        if(istatout(3,i).ne.izero) then
          nstatesrate=nstatesrate+ione
          itmp(nstatesrate,3)=i
        end if
      end do
      nstatesout=nstatestot+nstatesinc+nstatesrate
      do i=1,nstatestot
        itmp2(i)=itmp1(i,1)
      end do
      do i=1,nstatesinc
        itmp2(i+nstatestot)=itmp1(i,2)+nstatesmax
      end do
      do i=1,nstatesrate
        itmp2(i+nstatestot+nstatesinc)=itmp1(i,3)+2*nstatesmax
      end do
c
c...  create and open UCD file if UCD output is desired
c
      if(iucd.ne.izero) then
        call open_ucd(kucd,iprestress,nstep,ucdroot)
c******  can write header using itmp1 array
        call write_gauss_ucd_header(
c
c...  loop over element families
c
      do ifam=i,nvfamilies
        nelfamily=ivfamily(1,ifam)
        matmodel=ivfamily(2,ifam)
        indstate=ivfamily(3,ifam)
        nstate=infmatmod(2,matmodel)
        if(idout.gt.1) write(kw,700) ifam
c
        if(matmodel.eq.1) then
          call write_state_cmp(
     &     state(indstate),dstate(indstate),nelfamily,nstate,           ! elemfamily
     &     get_state_1,                                                 ! materl
     &     ngauss,                                                      ! eltype
     &     deltu,nstep,                                                 ! timdat
     &     itmp,                                                        ! ioopts
     &     idout,idsk,iucd,kw,kp,kucd)                                  ! ioinfo
c
c...  loop over number of state variables
c
      do i=1,nstatesmax
        if(istatout(1,i).ne.izero) then
          nout=izero
          do iel=1,numelt
            imat=infiel(2,iel)
            ietype=infiel(3,iel)
            ngauss=infetype(1,ietype)
            matmodel=infmat(1,imat)
            nstate=infmatmod(2,matmodel)
            indstate=infiel(5,iel)+matmodpt(i,matmodel)
            do l=1,ngauss
              indstateg=indstate+(l-1)*nstate
              if(ismatmod(i,matmodel).eq.izero) then
                call fill(stmp,zero,nstr)
              else
                call dcopy(nstr,state(1,indstateg),ione,stmp,ione)
              end if
              if(idout.gt.1) then
                nout=nout+1
                if(nout.eq.1.or.mod(nout,m50).eq.izero) then
                  write(kw,2000) i,(labels(j),j=1,nstr)
                  write(kw,*) ' '
                end if
                write(kw,3000) iel,l,(stmp(j),j=1,nstr)
              end if
              if(idsk.eq.ione) write(kp,1500) iel,l,(stmp(j),j=1,nstr)
              if(idsk.eq.itwo) write(kp) stmp
            end do
          end do
        end if
c
c...  output increments/rates, if desired.
c
        if(istatout(2,i).ne.izero) then
          tmult=one
          statedescr="   i n c r e m e n t "
          if(istatout(2,i).eq.ione) then
            if(delt.gt.zero) tmult=one/delt
            statedescr="   r a t e "
          end if
cdebug          write(6,*) "delt,nstep,tmult:",delt,nstep,tmult
          nout=izero
          do iel=1,numelt
            imat=infiel(2,iel)
            ietype=infiel(3,iel)
            ngauss=infetype(1,ietype)
            matmodel=infmat(1,imat)
            nstate=infmatmod(2,matmodel)
            indstate=infiel(5,iel)+matmodpt(i,matmodel)
            do l=1,ngauss
              indstateg=indstate+(l-1)*nstate
              call fill(stmp,zero,nstr)
              if(ismatmod(i,matmodel).ne.izero) call daxpy(nstr,tmult,
     &         dstate(1,indstateg),ione,stmp,ione)
              if(idout.gt.1) then
                nout=nout+1
                if(nout.eq.1.or.mod(nout,m50).eq.izero) then
                  write(kw,2500) statedescr,i,(labels(j),j=1,nstr)
                  write(kw,*) ' '
                end if
                write(kw,3000) iel,l,(stmp(j),j=1,nstr)
              end if
              if(idsk.eq.ione) write(kp,1500) iel,l,(stmp(j),j=1,nstr)
              if(idsk.eq.itwo) write(kp) stmp
            end do
          end do
        end if
      end do
      if(iucd.ne.izero) close(kucd)
      return
2000  format(///1x,' s t a t e   v a r i a b l e ',i3,///,
     &' elem #  gausspt',5x,5(a4,11x),a4)
2500  format(///1x,' s t a t e   v a r i a b l e ',a21,i3,///,
     &' elem #  gausspt',5x,5(a4,11x),a4)
3000  format(1x,i7,1x,i7,1x,6(1pe15.6))
1500  format(2i7,6e16.7)
      end
c
c version
c $Id: write_state_drv.f,v 1.1 2005/03/23 02:47:34 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
