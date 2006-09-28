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
      subroutine stiff_ss(
     & xl,iddmat,                                                       ! global
     & dmat,ien,iel,                                                    ! elemnt
     & gauss,sh,shj,nen,ngauss,nee,                                     ! eltype
     & s,                                                               ! stiff
     & getshape,bmatrix,                                                ! bbar
     & ierr,errstrng)                                                   ! errcode
c
c...computes the local stiffness matrix at the given integration points.
c   k=(b)t*d*b.
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
      integer iel,nen,ngauss,nee,ierr
      integer iddmat(nstr,nstr),ien(nen)
      character errstrng*(*)
      double precision xl(nsd,nen),dmat(nddmat,ngauss)
      double precision gauss(nsd+1,ngauss),sh(nsd+1,nen,ngauss)
      double precision shj(nsd+1,nen,ngauss),s(nee,nee)
c
c...  external routines
c
      external getshape,bmatrix
c
c...  local variables
c
      integer l
      double precision shd(4,nenmax,ngaussmax),b(6,3*nenmax)
      double precision shbar(4,nenmax),db(6,3*nenmax),det(ngaussmax)
      double precision dtmp(6,6)
      double precision vol
c
cdebug      integer idb,jdb
c
c...form shape functions for each integration point
c
cdebug      write(6,*) "Hello from stiff_ss_f!"
c
      call getshape(xl,sh,shj,shd,shbar,det,gauss,vol,iel,nen,ngauss,
     & ierr,errstrng)
      if(ierr.ne.izero) return
c
c...construct b matrix, then form intermediate db=dmat*b, and finally
c   the stiffness (b)t*dmat*b multiplied by appropriate weight for
c   integral over element
c
      do l=1,ngauss
        call bmatrix(b,shd(1,1,l),shbar,nen)
        call getmat(dtmp,dmat(1,l),iddmat,nstr,nddmat)
        call dsymm("l","l",nstr,nee,det(l),dtmp,nstr,b,nstr,zero,db,
     &    nstr)
c*        call dspmv("u",nstr,det(l),dmat(1,l),b,ione,zero,db,ione)
        call dgemm("t","n",nee,nee,nstr,one,b,nstr,db,nstr,one,s,nee)
cdebug        if(iel.eq.8) then
cdebug          write(6,*) "det:",det(l)
cdebug          do idb=1,nee
cdebug            write(6,*) "b:",(b(jdb,idb),jdb=1,nstr)
cdebug          end do
cdebug          write(6,*) "dmat:",(dmat(jdb,l),jdb=1,nddmat)
cdebug        end if
      end do
cdebug      do idb=1,nee
cdebug        write(6,*) "idb,s(idb,idb):",idb,s(idb,idb)
cdebug      end do
      return
      end
c
c version
c $Id: stiff_ss.f,v 1.6 2005/04/01 23:31:35 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
