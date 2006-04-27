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
      subroutine loadx(bconcforce,diforc,histry,idx,idhist,neq,
     & nconcflag,numnp,nhist,nstep,lastep,ierr,errstrng)
c
c...program to transfer differential forces into global load vector
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer neq,nconcflag,numnp,nhist,nstep,lastep,ierr
      integer idx(ndof,numnp),idhist(numnp)
      character errstrng*(*)
      double precision bconcforce(neq*nconcflag),diforc(ndof,numnp)
      double precision histry(nhist,lastep+1)
c
c...  local variables
c
      integer j,i,k,ihist
      double precision fac
c
cdebug      write(6,*) "Hello from loadx_f!"
c
      do j=1,numnp
        ihist=idhist(j)
        if((ihist.gt.nhist).or.(ihist.lt.izero)) then
          ierr=100
          errstrng="loadx"
          return
        end if
        do i=1,ndof
          k=idx(i,j)
          if(k.ne.izero) then
            fac=one
            if(ihist.gt.izero) fac=histry(ihist,nstep+1)
            bconcforce(k)=bconcforce(k)+diforc(i,j)*fac
          end if
        end do
      end do
      return
      end
c
c version
c $Id: loadx.f,v 1.3 2005/01/05 22:48:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
