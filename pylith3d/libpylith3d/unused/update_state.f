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
      integer infiel(7,numelt),infmat(3,numat),infmatmod(5,nmatmodmax)
      integer infetype(4,netypes)
      character errstrng*(*)
      double precision state(nstr,nstatesz),dstate(nstr,nstatesz)
c
c...  local variables
c
      integer matgpt,imat,matmodel,nmatel,nstate,ielg,iel
      integer ietype,indstate,indstateg,ngauss,l,nstr2,nstr3,nstr4
      double precision sub
c
cdebug      write(6,*) "Hello from update_state_f!"
c
      matgpt=ione
      sub=-one
      nstr2=itwo*nstr
      nstr3=ithree*nstr
      nstr4=ifour*nstr
c
c...  loop over material groups
c
      do imat=1,numat
        matmodel=infmat(1,imat)
        nmatel=infmat(2,imat)
        nstate=infmatmod(2,matmodel)
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
              call daxpy(nstr2,sub,state(1,indstateg),ione,
     &         dstate(1,indstateg),ione)
              call daxpy(nstr2,one,dstate(1,indstateg),ione,
     &         state(1,indstateg),ione)
              indstateg=indstateg+nstate
            end do
          end do
        else if(nstate.eq.ithree) then
          do ielg=matgpt,matgpt+nmatel-1
            iel=infiel(4,ielg)
            ietype=infiel(3,iel)
            indstate=infiel(5,iel)
            indstateg=indstate
            ngauss=infetype(1,ietype)
            do l=1,ngauss
              call daxpy(nstr2,sub,state(1,indstateg),ione,
     &         dstate(1,indstateg),ione)
              call daxpy(nstr3,one,dstate(1,indstateg),ione,
     &         state(1,indstateg),ione)
              indstateg=indstateg+nstate
            end do
          end do
        else if(nstate.eq.ifour) then
          do ielg=matgpt,matgpt+nmatel-1
            iel=infiel(4,ielg)
            ietype=infiel(3,iel)
            indstate=infiel(5,iel)
            indstateg=indstate
            ngauss=infetype(1,ietype)
            do l=1,ngauss
              call daxpy(nstr2,sub,state(1,indstateg),ione,
     &         dstate(1,indstateg),ione)
              call daxpy(nstr4,one,dstate(1,indstateg),ione,
     &         state(1,indstateg),ione)
              indstateg=indstateg+nstate
            end do
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
c $Id: update_state.f,v 1.1 2005/03/25 23:07:10 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
