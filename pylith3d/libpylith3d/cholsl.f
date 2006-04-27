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
      subroutine cholsl(a,n,np,p,b,x)
c
c...  routine to solve a square inverse problem given a
c     Cholesky-factored matrix a.
c     Taken from Numerical Recipes
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer n,np
      double precision a(np,np),b(n),p(n),x(n)
c
c...  local variables
c
      integer i,k
      double precision sum
c
cdebug      write(6,*) "Hello from cholsl_f!"
c
      do i=1,n
        sum=b(i)
        do k=i-1,1,-1
          sum=sum-a(i,k)*x(k)
        end do
        x(i)=sum/p(i)
      end do
      do i=n,1,-1
        sum=x(i)
        do k=i+1,n
          sum=sum-a(k,i)*x(k)
        end do
        x(i)=sum/p(i)
      end do
      return
      END
c
c version
c $Id: cholsl.f,v 1.2 2004/08/12 01:12:46 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
