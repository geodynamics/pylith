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
      subroutine update_state(state,dstate,infiel,infmat,infmatmod,
     & infetype,nstatesz,numelt,numat,ierr,errstrng)
c
c...program to update state variables after iteration convergence
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
      integer nstatesz,numelt,numat,ierr
      integer infiel(6,numelt),infmat(3,numat),infmatmod(5,nmatmodmax)
      integer infetype(4,netypes)
      character errstrng*(*)
      double precision state(nstr,nstatesz),dstate(nstr,nstatesz)
c
c...  local variables
c
      integer matgpt,imat,matmodel,nmatel,nstate,incstate,ielg,iel
      integer ietype,indstate,indstateg,ngauss,l
      double precision sub
c
cdebug      write(6,*) "Hello from update_state_f!"
c
      matgpt=ione
      sub=-one
c
c...  loop over material groups
c
      do imat=1,numat
        matmodel=infmat(1,imat)
        nmatel=infmat(2,imat)
        nstate=infmatmod(2,matmodel)
        incstate=nstr*nstate
c
c...  loop over elements in group, updating appropriate state variables
c
        if(nstate.eq.itwo) then
          do ielg=matgpt,matgpt+nmatel-1
            iel=infiel(4,ielg)
            ietype=infiel(3,iel)
            indstate=infiel(5,iel)
            indstateg=indstate
            ngauss=infetype(1,ietype)
            do l=1,ngauss
              call dcopy(nstr,dstate(1,indstateg),ione,
     &         state(1,indstateg),ione)
              call daxpy(nstr,sub,state(1,indstateg+nstr),ione,
     &         dstate(1,indstateg+nstr),ione)
              call daxpy(nstr,one,dstate(1,indstateg+nstr),ione,
     &         state(1,indstateg+nstr),ione)
            end do
            indstateg=indstateg+incstate
          end do
        else if(nstate.eq.ithree) then
          do ielg=matgpt,matgpt+nmatel-1
            iel=infiel(4,ielg)
            ietype=infiel(3,iel)
            indstate=infiel(5,iel)
            indstateg=indstate
            ngauss=infetype(1,ietype)
            do l=1,ngauss
              call dcopy(nstr,dstate(1,indstateg),ione,
     &         state(1,indstateg),ione)
              call daxpy(nstr,sub,state(1,indstateg+nstr),ione,
     &         dstate(1,indstateg+nstr),ione)
              call daxpy(nstr,one,dstate(1,indstateg+nstr),ione,
     &         state(1,indstateg+nstr),ione)
              call daxpy(nstr,one,dstate(1,indstateg+2*nstr),ione,
     &         state(1,indstateg+2*nstr),ione)
            end do
            indstateg=indstateg+incstate
          end do
        else if(nstate.eq.ifour) then
          do ielg=matgpt,matgpt+nmatel-1
            iel=infiel(4,ielg)
            ietype=infiel(3,iel)
            indstate=infiel(5,iel)
            indstateg=indstate
            ngauss=infetype(1,ietype)
            do l=1,ngauss
              call dcopy(nstr,dstate(1,indstateg),ione,
     &         state(1,indstateg),ione)
              call daxpy(nstr,sub,state(1,indstateg+nstr),ione,
     &         dstate(1,indstateg+nstr),ione)
              call daxpy(nstr,one,dstate(1,indstateg+nstr),ione,
     &         state(1,indstateg+nstr),ione)
              call daxpy(nstr,one,dstate(1,indstateg+2*nstr),ione,
     &         state(1,indstateg+2*nstr),ione)
              call daxpy(nstr,one,dstate(1,indstateg+3*nstr),ione,
     &         state(1,indstateg+3*nstr),ione)
            end do
            indstateg=indstateg+incstate
          end do
        else
          ierr=102
          errstrng="update_state"
          return
        end if
        matgpt=matgpt+nmatel
      end do
      return
      end
c
c version
c $Id: update_state.f,v 1.4 2004/07/12 21:07:20 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
