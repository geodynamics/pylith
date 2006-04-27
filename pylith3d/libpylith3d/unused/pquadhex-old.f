c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine pquadhex(sh,shj,gauss,ngauss,nen,nsd,nenmax,ngaussmax,
     & intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a quadratic (20-node)
c    hexahedron.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nsd,nenmax,ngaussmax,intord
      integer ngauss,nen
      double precision sh(nsd+1,nenmax,ngaussmax)
      double precision shj(nsd+1,nenmax,ngaussmax)
      double precision gauss(nsd+1,ngaussmax)
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
      double precision r(20),s(20),t(20)
      data r/-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0,
     &        0d0, 1d0, 0d0,-1d0, 0d0, 1d0, 0d0,-1d0,
     &       -1d0, 1d0, 1d0,-1d0/
      data s/-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,
     &       -1d0, 0d0, 1d0, 0d0,-1d0, 0d0, 1d0, 0d0,
     &       -1d0,-1d0, 1d0, 1d0/
      data t/ 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0,
     &        1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0,
     &        0d0, 0d0, 0d0, 0d0/
c
c...  intrinsic functions
c
      intrinsic abs,sqrt,dble
c
c...  user-defined functions
c
      double precision gquad,dgquad
c
c...  local variables
c
      integer i,l,l1,l2,l3,nshsize,ngssize
      double precision g1,w1,w2,rr,ss,tt,rrw,ssw,ttw,drr,dss,dtt,beta
      double precision betai
c
c...  function definitions
c
      gquad(beta,betai)=(half*(one-abs(betai))+half)*(one+beta*betai+
     & (abs(betai)-one)*beta*beta)
      dgquad(beta,betai)=half*betai-two*beta*(one-abs(betai))
c
c...  definitions
c
      nshsize=(nsd+1)*nenmax*ngaussmax
      ngssize=(nsd+1)*ngaussmax
c
c...  Quadratic hex definition
c
      nen=20
      ngauss=ieight
      do l=1,ngauss
        gauss(1,l)=r(l)*root3i
        gauss(2,l)=s(l)*root3i
        gauss(3,l)=t(l)*root3i
        gauss(4,l)=one
      end do
      if(intord.ne.2) then
        ngauss=27
        g1=sqrt(three/five)
        w1=five/nine
        w2=eight/nine
        l=0
        do l3=1,3
          tt=dble(l3-2)
          ttw=tt*w1+(abs(tt)-one)*w2
          do l2=1,3
            ss=dble(l2-2)
            ssw=ss*w1+(abs(ss)-one)*w2
            do l1=1,3
              rr=dble(l1-2)
              rrw=rr*w1+(abs(rr)-one)*w2
              l=l+1
              gauss(1,l)=rr*g1
              gauss(2,l)=ss*g1
              gauss(3,l)=tt*g1
              gauss(4,l)=rrw*ssw*ttw
            end do
          end do
        end do
      end if
      do l=1,ngauss
        do i=1,nen
          rr=gquad(gauss(1,l),r(i))
          ss=gquad(gauss(2,l),s(i))
          tt=gquad(gauss(3,l),t(i))
          drr=dgquad(gauss(1,l),r(i))
          dss=dgquad(gauss(2,l),s(i))
          dtt=dgquad(gauss(3,l),t(i))
          sh(4,i,l)=rr*ss*tt
          sh(1,i,l)=drr*ss*tt
          sh(2,i,l)=dss*rr*tt
          sh(3,i,l)=dtt*rr*ss
        end do
        do i=1,4
          sh(i,1,l)=sh(i,1,l)-half*(sh(i,9,l)+sh(i,12,l)+ sh(i,17,l))
          sh(i,2,l)=sh(i,2,l)-half*(sh(i,9,l)+sh(i,10,l)+sh(i,18,l))
          sh(i,3,l)=sh(i,3,l)-half*(sh(i,10,l)+sh(i,11,l)+sh(i,19,l))
          sh(i,4,l)=sh(i,4,l)-half*(sh(i,11,l)+sh(i,12,l)+sh(i,20,l))
          sh(i,5,l)=sh(i,5,l)-half*(sh(i,13,l)+sh(i,16,l)+sh(i,17,l))
          sh(i,6,l)=sh(i,6,l)-half*(sh(i,13,l)+sh(i,14,l)+sh(i,18,l))
          sh(i,7,l)=sh(i,7,l)-half*(sh(i,14,l)+sh(i,15,l)+sh(i,19,l))
          sh(i,8,l)=sh(i,8,l)-half*(sh(i,15,l)+sh(i,16,l)+sh(i,20,l))
        end do
      end do
      call dcopy(nshsize,sh,ione,shj,ione)
c
      return
      end
c
c version
c $Id: pquadhex-old.f,v 1.1 2004/07/06 19:20:28 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
