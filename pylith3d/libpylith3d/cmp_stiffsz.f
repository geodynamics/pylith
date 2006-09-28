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
      subroutine cmp_stiffsz(neq,lm,lmx,numelv,iwork,numsn,nen,
     & ierr,errstrng)
c
c      program to compute an upper bound on the number of nonzero
c      entries in the stiffness matrix.  Duplicate contributions from
c      different elements are not considered, so this is quite an
c      overestimate.
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
      integer neq,numelv,iwork,numsn,nen,ierr
      integer lm(ndof*nen,numelv),lmx(ndof*nen,numelv)
      character errstrng*(*)
c
c...  intrinsic functions
c
      intrinsic abs
c
c...  local variables
c
      integer iel,nee,neeq,i,irow,j,icol
c
cdebug      write(6,*) "Hello from cmp_stiffsz_f!"
c
      iwork=izero
      nee=ndof*nen
      neeq=nee
      if(numsn.ne.izero) neeq=itwo*nee
      do iel=1,numelv
c
        do i=1,neeq
          if(i.le.nee) then
            irow=lm(i,iel)
          else
            irow=abs(lmx(i-nee,iel))
          end if
          if(irow.ne.izero) then
            do j=1,neeq
              if(j.le.nee) then
                icol=lm(j,iel)
              else
                icol=abs(lmx(j-nee,iel))
              end if
              if(icol.ne.izero.and.icol.ne.irow) iwork=iwork+ione
            end do
          end if
        end do
      end do
      iwork=iwork+neeq*neeq
      return
      end
c
c version
c $Id: cmp_stiffsz.f,v 1.5 2005/04/16 00:39:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
