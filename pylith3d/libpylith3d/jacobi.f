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
      subroutine jacobi(a,n,np,d,v,nrot)
c
c...routine to compute the eigenvalues and eigenvectors of matrix a,
c   using Jacobi rotations.  Taken from Numerical Recipes.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "rconsts.inc"
      integer nmax
      parameter (nmax=100)
      double precision eps,tval,gval
      parameter(eps=1.0d-18,tval=0.2d0,gval=100.0d0)
c
c...  subroutine arguments
c
      integer n,np,nrot
      double precision a(np,np),d(np),v(np,np)
c
c...  intrinsic functions
c
      intrinsic abs,sqrt,dble
c
c...  local variables
c
      integer i,ip,iq,j
      double precision c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)
c
cdebug      write(6,*) "Hello from jacobi_f!"
c
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=zero
11      continue
        v(ip,ip)=one
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=zero
13    continue
      nrot=0
      do 24 i=1,50
        sm=zero
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(abs(sm).lt.eps)return
        if(i.lt.4)then
          tresh=tval*sm/(dble(n)*dble(n))
        else
          tresh=zero
        end if
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=gval*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     *         .and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=zero
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=half*h/a(ip,iq)
                t=one/(abs(theta)+sqrt(one+theta**2))
                if(theta.lt.zero)t=-t
              end if
              c=one/sqrt(one+t*t)
              s=t*c
              tau=s/(one+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=zero
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            end if
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=zero
23      continue
24    continue
      pause '50 iterations should never happen'
      return
      end
c
c version
c $Id: jacobi.f,v 1.3 2004/08/12 01:33:20 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
