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
      subroutine elas_mat_cmp(
     & dmat,infiel,iddmat,ndmatsz,numelt,                               ! elemnt
     & prop,infmat,nprop,matgpt,elas_mat,                               ! materl
     & infetype,                                                        ! eltype
     & ierr,errstrng)                                                   ! errcode
c
c...  compute subroutine to form the d-matrix for the elastic solution
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nshape.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer ndmatsz,numelt,nprop,matgpt,ierr
      integer infmat(6),infiel(6,numelt),infetype(4,netypes)
      integer iddmat(nstr,nstr)
      character errstrng*(*)
      double precision dmat(nddmat,ndmatsz),prop(nprop)
c
c...  external routines
c
      external elas_mat
c
c...  local variables
c
      integer nmatel,imatvar,iel,inddmat0,ietype,ngauss,ng,inddmatg
      integer l,ind,inddmat
c
cdebug      write(6,*) "Hello from elas_mat_cmp_f!"
c
      nmatel=infmat(2)
      imatvar=infmat(4)
c
c...  compute d-matrix for first element in the group.  If imatvar is
c     zero, there is only one d-matrix for the group, and that is all
c     that needs to be done.  Otherwise, since this is the elastic
c     solution, the same matrix can be copied into each slot.
c
      iel=infiel(4,matgpt)
      inddmat0=infiel(6,iel)
      ietype=infiel(3,iel)
      call elas_mat(dmat(1,inddmat0),prop,iddmat,nprop,ierr,errstrng)
      if(ierr.ne.izero) return
      ngauss=infetype(1,ietype)
      ng=ngauss
      if(imatvar.eq.izero) ng=ngaussmax
      inddmatg=inddmat0
      do l=2,ng
        inddmatg=inddmatg+nddmat
        call dcopy(nddmat,dmat(1,inddmat0),ione,dmat(1,inddmatg),ione)
      end do
c
c...  loop over elements in group if there is material property
c     variation for the material type.
c
      if(imatvar.ne.0) then
        do ind=matgpt+1,matgpt+nmatel-1
          iel=infiel(4,ind)
          ietype=infiel(3,iel)
          inddmat=infiel(6,iel)
          ngauss=infetype(1,ietype)
          inddmatg=inddmat
          do l=1,ngauss
            inddmatg=inddmat+nddmat
            call dcopy(nddmat,dmat(1,inddmat0),ione,dmat(1,inddmatg),
     &       ione)
          end do
        end do
      end if
      return
      end
c
c version
c $Id: elas_mat_cmp.f,v 1.1 2004/06/25 21:41:22 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
