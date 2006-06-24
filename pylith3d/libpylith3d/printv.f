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
      subroutine printv(a,b,id,idx,neq,numnp,iout,idout,kw)
c
c...program to print one or two vectors of length neq, along with
c   data on correlation between equation number and node
c
c     if iout = 0, print array d after factorzation a=(u)t*d*u
c        iout = 1, print global load vectors for debugging
c        iout = 2, print diagonals of stiffness matrix
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
c
c...  subroutine arguments
c
      integer neq,numnp,iout,idout,kw
      integer id(ndof,numnp),idx(ndof,numnp)
      double precision a(neq),b(neq)
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer npage,n,i,j,iflag,idof,node
c
cdebug      write(6,*) "Hello from printv_f!"
c
      if(idout.lt.2) return
      npage=50
      do n=1,neq
c
c...find node and dof associated with equation number
c
        do i=1,ndof
          do j=1,numnp
            iflag=0
            if(id(i,j).eq.n) goto 50
            if(idx(i,j).eq.n) then
              iflag=1
              goto 50
            end if
          end do
        end do
   50   idof=i
        node=j
c
c...print out the appropriate vectors
c
        if(iout.eq.0) then
          if((n.eq.1).or.(mod(n,npage).eq.0)) write(kw,1500)
          if(iflag.eq.0) write(kw,2000) n,node,idof,a(n)
          if(iflag.eq.1) write(kw,3000) n,node,idof,a(n)
        else if(iout.eq.1) then
          if((n.eq.1).or.(mod(n,npage).eq.0)) write(kw,1000)
          if(iflag.eq.0) write(kw,2000) n,node,idof,a(n),b(n),
     &     b(n)-a(n)
          if(iflag.eq.1) write(kw,3000) n,node,idof,a(n),b(n),
     &     b(n)-a(n)
        else if(iout.eq.2) then
          if((n.eq.1).or.(mod(n,npage).eq.0)) write(kw,1700)
          if(iflag.eq.0) write(kw,2000) n,node,idof,a(n)
          if(iflag.eq.1) write(kw,3000) n,node,idof,a(n)
        end if
      end do
      return
 1000 format(//' global force vector follows'//1x,'  eq# ',' node',
     & '    dof  ',5x,' current b',7x,'applied b',5x,'difference'/)
 1500 format(//' diagonal elements of stiffness matrix:'///
     & 2x,' neq ',' node ','  dof ',7x, '    diag'//)
 1700 format(//' diagonal elements of stiffness matrix:'///
     & 2x,' neq ',' node ','  dof ',7x, '    diag'//)
 2000 format(1x,3(i7,1x ),1x,3(1x,1pe15.8))
 3000 format(1x,3(i7,'x'),1x,3(1x,1pe15.8))
      end
c
c version
c $Id: printv.f,v 1.2 2004/07/05 19:42:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
