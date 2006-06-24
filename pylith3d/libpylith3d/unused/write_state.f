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
      subroutine write_state(
     & state,dstate,infiel,nstatesz,numelt,                             ! elemnt
     & infmat,infmatmod,ismatmod,numat,                                 ! materl
     & infetype,                                                        ! eltype
     & delt,nstep,                                                      ! timdat
     & istatout,                                                        ! ioopts
     & idout,idsk,kw,kp)                                                ! ioinfo
c
c...  program to print stress and strain
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
      integer nstatesz,numelt,numat,nstep,idout,idsk,kw,kp
      integer infiel(7,numelt),infmat(3,numat),infmatmod(5,nmatmodmax)
      integer ismatmod(nstatesmax,nmatmodmax),infetype(4,netypes)
      integer istatout(2,nstatesmax)
      double precision delt
      double precision state(nstr,nstatesz),dstate(nstr,nstatesz)
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
      m50=50
c
c...  define pointer array for location of state variables for each
c     material model.
c
      do i=1,nmatmodmax
        matmodpt(1,i)=izero
        do j=2,nstatesmax
          matmodpt(j,i)=matmodpt(j-1,i)+ismatmod(j-1,i)
        end do
      end do
cdebug      do idb=1,nstatesz
cdebug        write(6,*) "state:",(state(jdb,idb),jdb=1,nstr)
cdebug      end do
cdebug      do idb=1,nstatesz
cdebug        write(6,*) "dstate:",(dstate(jdb,idb),jdb=1,nstr)
cdebug      end do
cdebug      do idb=1,nmatmodmax
cdebug        write(6,*) "matmodpt:",(matmodpt(jdb,idb),jdb=1,nstatesmax)
cdebug      end do
cdebug      do idb=1,nmatmodmax
cdebug        write(6,*) "ismatmod:",(ismatmod(jdb,idb),jdb=1,nstatesmax)
cdebug      end do
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
c $Id: write_state.f,v 1.1 2005/03/25 22:52:39 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
