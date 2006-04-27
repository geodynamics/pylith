c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams
c  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine write_state_cmp(
     & state,dstate,nelfamily,nstate,ielg,                              ! elemfamily
     & get_state,                                                       ! materl
     & ngauss,                                                          ! eltype
     & delti,nstep,                                                     ! timdat
     & istatout,istatoutc,nstatout,nstatestot,nout,npts,                ! ioopts
     & statemin,statemax,ibyteg,                                        ! binucd
     & idout,idsk,iucd,kw,kp,kucd)                                      ! ioinfo
c
c...  computation routine to update state variables within an element family
c
#include "binio.h"
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "materials.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nelfamily,nstate,ielg,ngauss,nstep,nstatestot,nout,npts
      integer idout,idsk,iucd,kw,kp,kucd
      integer istatout(nstatesmax,3),istatoutc(nstatestot),nstatout(3)
      integer ibyteg(3*nstatesmax)
      double precision state(nstate,ngauss,nelfamily)
      double precision dstate(nstate,ngauss,nelfamily)
      double precision statemin(3*nstatesmax),statemax(3*nstatesmax)
      double precision delti
c
c...  external routines
c
      external get_state
c
c...  included dimension and type statements
c
      include "labels_dim.inc"
c
c...  local constants
c
      integer npage
      data npage/50/
c
c...  local variables
c
      integer ielga,ielgp,ielf,j,l
      integer floatlen,istate
      double precision sout(3*nstatesmax)
c
cdebug      integer idb
c
c...  included variable definitions
c
      include "labels_def.inc"
c
cdebug      write(6,*) "Hello from write_state_cmp_f!"
c
      ielga=ielg
      ielgp=ielg
      floatlen=ifour
      call fill(sout,zero,3*nstatesmax)
c
c...  output to human-readable ascii file
c
      if(idout.gt.1) then
        do ielf=1,nelfamily
          do l=1,ngauss
            nout=nout+ione
            call get_state(state(1,l,ielf),dstate(1,l,ielf),sout,nstate)
            call daxpy(nstatesmax,delti,sout(nstatesmax+ione),ione,
     &       sout(2*nstatesmax+ione),ione)
cdebug            write(6,*) "sout:",(sout(idb),idb=1,3*nstatesmax)
            if(nout.eq.1.or.mod(nout,npage).eq.izero) then
              write(kw,700) (labels(istatoutc(j)),j=1,nstatestot)
              write(kw,"(/)")
            end if
            write(kw,710) ielga,l,(sout(istatoutc(j)),j=1,nstatestot)
          end do
          ielga=ielga+ione
        end do
      end if
c
c...  output to ascii plot file
c
      if(idsk.eq.1) then
        do ielf=1,nelfamily
          do l=1,ngauss
            call get_state(state(1,l,ielf),dstate(1,l,ielf),sout,nstate)
            call daxpy(nstatesmax,delti,sout(nstatesmax+ione),ione,
     &       sout(2*nstatesmax+ione),ione)
            write(kw,720) ielgp,l,(sout(istatoutc(j)),j=1,nstatestot)
          end do
          ielgp=ielgp+ione
        end do
      end if
c
c...  output to binary plot file
c
      if(idsk.eq.2) then
        do ielf=1,nelfamily
          do l=1,ngauss
            call get_state(state(1,l,ielf),dstate(1,l,ielf),sout,nstate)
            call daxpy(nstatesmax,delti,sout(nstatesmax+ione),ione,
     &       sout(2*nstatesmax+ione),ione)
            write(kp) (sout(istatoutc(j)),j=1,nstatestot)
          end do
        end do
      end if
c
c...  output to UCD file
c
      if(iucd.eq.1) then
        do ielf=1,nelfamily
          do l=1,ngauss
            npts=npts+ione
            call get_state(state(1,l,ielf),dstate(1,l,ielf),sout,nstate)
            call daxpy(nstatesmax,delti,sout(nstatesmax+ione),ione,
     &       sout(2*nstatesmax+ione),ione)
cdebug            write(6,*) "sout:",(sout(idb),idb=1,3*nstatesmax)
            write(kucd,720) npts,(sout(istatoutc(j)),j=1,nstatestot)
          end do
        end do
      else if(iucd.eq.2) then
        do ielf=1,nelfamily
          do l=1,ngauss
            call get_state(state(1,l,ielf),dstate(1,l,ielf),sout,nstate)
            call daxpy(nstatesmax,delti,sout(nstatesmax+ione),ione,
     &       sout(2*nstatesmax+ione),ione)
            do j=1,nstatestot
              istate=istatoutc(j)
              write(kucd,BINARY_POS=ibyteg(j)) real(sout(istate))
              statemin(istate)=min(statemin(istate),sout(istate))
              statemax(istate)=max(statemax(istate),sout(istate))
              ibyteg(j)=ibyteg(j)+floatlen
            end do
          end do
        end do
      end if
c
 700  format(//,2x,"elmt",4x,"gspt",100(:5x,a11))
 710  format(2(1x,i7),100(:2x,1pe15.8))
 720  format(i7,100(:2x,1pe15.8))
      return
      end
c
c version
c $Id: write_state_cmp.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c End of file 
