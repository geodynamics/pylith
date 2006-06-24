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
      subroutine skcomp(x,d,skew,idslp,ipslp,ipstrs,numsn,numnp,nstep,
     & lgdefp,ierr,errstrng)
c
c...  subroutine to compute skew angles for a given load vector and
c     fault orientation
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
      double precision eps,rdeg
      parameter(eps=1.0d-8,rdeg=180.0d0/pi)
c
c...  subroutine arguments
c
      integer ipstrs,numsn,numnp,nstep,lgdefp,ierr
      integer idslp(numsn),ipslp(npdim,numsn)
      double precision x(nsd,numnp),d(ndof,numnp),skew(nskdim,numnp)
      character errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic acos,sign,dble,sqrt,min,sin,cos,asin,tan,atan
c
c...  local variables
c
      integer ma,nd,ncvm,ldtmp,nn,n,i
      double precision rm,sdmin,xd,yd,zd,amag,bmag,sgn,sk1
      double precision dnm,rat
      double precision a(3),xx(5),yy(5),sd(5),xl(3,5),cov(3,3)
c
cdebug      write(6,*) "Hello from skcomp_f!"
c
      ma=ithree
      nd=ifive
      if(numsn.eq.ifour) nd=ifour
      if(numsn.eq.ithree) nd=ithree
      if(numsn.lt.ithree) then
        ierr=112
        errstrng="skcomp"
        return
      end if
      ncvm=ithree
      ldtmp=lgdefp
      if(ipstrs.eq.ione.and.nstep.eq.izero) ldtmp=izero
      rm=one
      if(ldtmp.eq.izero) rm=zero
c
c...  loop over slippery nodes, and compute normal for each plane
c
      do nn=1,numsn
        n=idslp(nn)
        sdmin=big
c
c...  localize coordinates and set up values for least-squares fit to
c     a plane
c
        xl(1,nd)=x(1,n)+rm*d(1,n)
        xl(2,nd)=x(2,n)+rm*d(2,n)
        xl(3,nd)=x(3,n)+rm*d(3,n)
        xx(nd)=dble(nd)
        yy(nd)=one
        do i=1,nd-1
          xl(1,i)=x(1,ipslp(i,nn))+rm*d(1,ipslp(i,nn))
          xl(2,i)=x(2,ipslp(i,nn))+rm*d(2,ipslp(i,nn))
          xl(3,i)=x(3,ipslp(i,nn))+rm*d(3,ipslp(i,nn))
          xd=xl(1,i)-xl(1,nd)
          yd=xl(2,i)-xl(2,nd)
          zd=xl(3,i)-xl(3,nd)
          xx(i)=dble(i)
          sd(i)=sqrt(xd*xd+yd*yd+zd*zd)
          yy(i)=one
          sdmin=min(sdmin,sd(i))
        end do
        sd(nd)=half*sdmin
        call lfit(xx,yy,sd,nd,a,ma,cov,ncvm,xl,ierr,errstrng)
        if(ierr.ne.izero) return
        amag=one/sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        a(1)=a(1)*amag
        a(2)=a(2)*amag
        a(3)=a(3)*amag
        amag=one/sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        bmag=zero
        sgn=one
        if(a(3).lt.zero) sgn=-one
        do i=1,ndof
          a(i)=sgn*a(i)*amag
          if(abs(a(i)).lt.eps) a(i)=zero
          bmag=bmag+a(i)*a(i)
        end do
        bmag=one/sqrt(bmag)
        do i=1,ndof
          a(i)=a(i)*bmag
          sgn=sign(one,a(i))
          if(abs(a(i)).gt.one) a(i)=sgn
        end do
c
c...  new z-axis is perpendicular to fault plane
c
        sk1=acos(a(3))
        if(abs(a(3)-one).lt.eps) then
          skew(1,n)=zero
          skew(2,n)=zero
        else if(abs(a(1)).gt.eps) then
          skew(1,n)=atan(a(2)/a(1))
          dnm=sin(skew(1,n))+cos(skew(1,n))
          rat=(a(1)+a(2))/dnm
          sgn=sign(one,rat)
          if(abs(rat).gt.one) rat=sgn
          skew(2,n)=asin(rat)
        else if(abs(a(3)).gt.eps) then
          skew(2,n)=sk1
          dnm=a(3)*tan(skew(2,n))
          rat=a(2)/dnm
          sgn=sign(one,rat)
          if(abs(rat).gt.one) rat=sgn
          skew(1,n)=asin(rat)
        else
          skew(2,n)=half*pi
          skew(1,n)=asin(a(2))
        end if
c*        write(17,710) n,a(1),a(2),a(3),rdeg*skew(1,n),rdeg*skew(2,n)
      end do
c*      write(17,720)
c*      call flush(17)
710   format(i6,5(2x,1pe12.5))
720   format(///)
      return
      end
c
c version
c $Id: skcomp.f,v 1.4 2004/08/12 02:26:24 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
