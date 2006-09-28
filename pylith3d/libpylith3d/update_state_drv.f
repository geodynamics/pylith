c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
      integer ifam,nelfamily,matmodel,indstate,nstate
c
cdebug      write(6,*) "Hello from update_state_drv_f!"
c
c
c...  loop over element families
c
      do ifam=1,nvfamilies
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
c $Id: update_state_drv.f,v 1.3 2005/04/05 23:02:40 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
