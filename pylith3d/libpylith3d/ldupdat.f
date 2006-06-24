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
      subroutine ldupdat(d,dx,tfault,dl,xl,ien,lmx,lmf,
     & nen,numnp,numfn,numslp)
c
c...subroutine to update local coordinates for large deformation
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
      integer nen,numnp,numfn,numslp
      integer ien(nen),lmx(ndof,nen),lmf(nen)
      double precision d(ndof,numnp),dx(ndof,numnp),tfault(ndof,numfn)
      double precision dl(ndof*nen),xl(nsd*nen)
c
c...  local variables
c
c
cdebug      write(6,*) "Hello from ldupdat_f"
c
      call ldisp(dl,d,ien,nen,numnp)
      if(numfn.ne.0) call adfldp(dl,lmf,tfault,nen,numfn)
      if(numslp.ne.0) call addsn(dl,dx,ien,lmx,nen,numnp)
      call daxpy(nen*nsd,one,dl,ione,xl,ione)
clater      if(iopt.eq.1.and.ldtmp.ne.0) then
clater        if(nsd.eq.ndof) then
clater          call daxpy(nen*nsd,one,dl,ione,xl,ione)
c
c...  Ugly (temporary) kludge to handle case where ndof.ne.nsd.
c     I could use this same setup for all cases, but it's probably
c     slower, and the special case only applies to out-of-plane
c     problems.
c
clater        else
clater          do i=1,nsd
clater            call daxpy(nen,one,dl(i),ndof,xl(i),nsd)
clater          end do
clater        end if
clater      end if
      return
      end
c
c version
c $Id: ldupdat.f,v 1.3 2004/08/12 01:34:40 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
