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
      subroutine update_state_drv(state,dstate,ivfamily,infmatmod,
     & nvfamilies,numelv,nstatesz,ngauss,ierr,errstrng)
c
c...program to update state variables after iteration convergence.
c   After updating, the values in state should represent current total
c   values, and those in dstate should represent changes since the
c   previous time step.
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
      integer nvfamilies,numelv,nstatesz,ngauss,ierr
      integer ivfamily(5,nvfamilies),infmatmod(6,nmatmodmax)
      character errstrng*(*)
      double precision state(nstatesz),dstate(nstatesz)
c
c...  external routines
c
      include "update_ext.inc"
c
c...  local variables
c
      integer ielg,ifam,nelfamily,matmodel,indstate,nstate
c
cdebug      write(6,*) "Hello from update_state_f!"
c
      ielg=ione
c
c...  loop over element families
c
      do ifam=1,nelfamily
        nelfamily=ivfamily(1,ifam)
        matmodel=ivfamily(2,ifam)
        indstate=ivfamily(3,ifam)
        nstate=infmatmod(2,matmodel)
        if(matmodel.eq.1) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_1)
        else if(matmodel.eq.2) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_2)
        else if(matmodel.eq.3) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_3)
        else if(matmodel.eq.4) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_4)
        else if(matmodel.eq.5) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_5)
        else if(matmodel.eq.6) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_6)
        else if(matmodel.eq.7) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_7)
        else if(matmodel.eq.8) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_8)
        else if(matmodel.eq.9) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_9)
        else if(matmodel.eq.10) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_10)
        else if(matmodel.eq.11) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_11)
        else if(matmodel.eq.12) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_12)
        else if(matmodel.eq.13) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_13)
        else if(matmodel.eq.14) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_14)
        else if(matmodel.eq.15) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_15)
        else if(matmodel.eq.16) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_16)
        else if(matmodel.eq.17) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_17)
        else if(matmodel.eq.18) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_18)
        else if(matmodel.eq.19) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_19)
        else if(matmodel.eq.20) then
          call update_state_cmp(state(indstate),dstate(indstate),
     &     nelfamily,nstate,ngauss,update_state_20)
        else
          ierr=101
          errstrng="update_state_drv"
          return
        end if
      end do
      return
      end
c
c version
c $Id: update_state_drv.f,v 1.1 2005/03/22 22:23:26 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
