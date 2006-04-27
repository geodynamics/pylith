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
      subroutine pquadtet(sh,shj,gauss,infetype,intord)
c
c... Subroutine to compute shape functions in natural coordinates,
c    integration points, and weights for a quadratic tetrahedron.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer intord
      integer infetype(4)
      double precision sh(nsd+1,nenmax,ngaussmax)
      double precision shj(nsd+1,nenmax,ngaussmax)
      double precision gauss(nsd+1,ngaussmax)
c
c...  local constants
c
c*      double precision r(10),s(10),t(10),u(10)
c*      data r/ 1d0, 0d0, 0d0, 0d0,0.5d0,0.0d0,0.5d0,0.5d0,0.0d0,0.0d0/
c*      data s/ 0d0, 1d0, 0d0, 0d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0,0.0d0/
c*      data t/ 0d0, 0d0, 1d0, 0d0,0.0d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0/
c*      data u/ 0d0, 0d0, 0d0, 1d0,0.0d0,0.0d0,0.0d0,0.5d0,0.5d0,0.5d0/
c
c...  intrinsic functions
c
      intrinsic sqrt,max
c
c...  user-defined functions
c
c
c...  local variables
c
c*      integer nen,ngauss,nec,nee,i,l,nshsize,ngssize
      integer nen,ngauss,nec,nee,l,nshsize,ngssize
c*      double precision g1,g2,rr,ss,tt,uu,drr,dss,dtt,duu
      double precision g1,g2,rr,ss,tt,uu
      double precision tetvol
cdebug      integer idb,jdb
cdebug      double precision dbsum
c
cdebug      write(6,*) "Hello from pquadtet_f!"
c
c...  definitions
c
      tetvol=sixth
      nshsize=(nsd+1)*nenmax*ngaussmax
      ngssize=(nsd+1)*ngaussmax
c
c...  Quadratic tet definition
c
      nen=10
      nec=nsd*nen
      nee=ndof*nen
      if(intord.eq.2) then
        ngauss=ione
        gauss(1,1)=fourth
        gauss(2,1)=fourth
        gauss(3,1)=fourth
        gauss(4,1)=tetvol
      else
        ngauss=ifour
        g1=(five-sqrt(five))/20.0d0
        g2=(five+three*sqrt(five))/20.0d0
        do l=1,ngauss
          gauss(1,l)=g1
          gauss(2,l)=g1
          gauss(3,l)=g1
          if(l.ne.ngauss) gauss(l,l)=g2
          gauss(4,l)=fourth*tetvol
        end do
      end if
c
      infetype(1)=ngauss
      infetype(2)=nen
      infetype(3)=nec
      infetype(4)=nee
c
      do l=1,ngauss
        rr=gauss(1,l)
        ss=gauss(2,l)
        tt=gauss(3,l)
        uu=g1
        if(l.eq.4) uu=g2
        if(intord.eq.itwo) uu=fourth
cdebug        write(6,*) "l:",l
cdebug        write(6,*) "rr,ss,tt,uu:",rr,ss,tt,uu
        sh(4,1,l)=rr*(two*rr-one)
        sh(1,1,l)=four*rr-one
        sh(2,1,l)=zero
        sh(3,1,l)=zero
        sh(4,2,l)=ss*(two*ss-one)
        sh(1,2,l)=zero
        sh(2,2,l)=four*ss-one
        sh(3,2,l)=zero
        sh(4,3,l)=tt*(two*tt-one)
        sh(1,3,l)=zero
        sh(2,3,l)=zero
        sh(3,3,l)=four*tt-one
        sh(4,4,l)=uu*(two*uu-one)
        sh(1,4,l)=-(four*uu-one)
        sh(2,4,l)=-(four*uu-one)
        sh(3,4,l)=-(four*uu-one)
        sh(4,5,l)=four*rr*ss
        sh(1,5,l)=four*ss
        sh(2,5,l)=four*rr
        sh(3,5,l)=zero
        sh(4,6,l)=four*ss*tt
        sh(1,6,l)=zero
        sh(2,6,l)=four*tt
        sh(3,6,l)=four*ss
        sh(4,7,l)=four*rr*tt
        sh(1,7,l)=four*tt
        sh(2,7,l)=zero
        sh(3,7,l)=four*rr
        sh(4,8,l)=four*rr*uu
        sh(1,8,l)=four*(uu-rr)
        sh(2,8,l)=-(four*rr)
        sh(3,8,l)=-(four*rr)
        sh(4,9,l)=four*ss*uu
        sh(1,9,l)=-(four*ss)
        sh(2,9,l)=four*(uu-ss)
        sh(3,9,l)=-(four*ss)
        sh(4,10,l)=four*tt*uu
        sh(1,10,l)=-(four*tt)
        sh(2,10,l)=-(four*tt)
        sh(3,10,l)=four*(uu-tt)
cdebug        dbsum=zero
cdebug        do idb=1,nen
cdebug          write(6,*) "sh:",(sh(jdb,idb,l),jdb=1,4)
cdebug          dbsum=dbsum+sh(4,idb,l)
cdebug        end do
cdebug        write(6,*) "l, sum sh:",l,dbsum
cdebug        write(6,*) "gauss:",(gauss(idb,l),idb=1,4)
      end do
      call dcopy(nshsize,sh,ione,shj,ione)
c
      return
      end
c
c version
c $Id: pquadtet-old.f,v 1.1 2004/08/12 14:39:02 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
