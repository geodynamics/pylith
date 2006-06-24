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
      subroutine getjac(x,xs,det,shj,nen,iel,ierr,errstrng)
c
c...  subroutine to compute the jacobian determinant given the element
c     coordinates and the shape functions in natural coordinates.
c
c       shj(1,nen),sh(2,nen),sh(3,nen) = x,y,and z derivatives
c                                        of shape functions
c       xs(nsd,nsd)                   = jacobian matrix
c       det                           = determinant of jacobian matrix
c       x(nsd,nen)                    = local nodal coordinates
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
      integer nen,iel,ierr
      character errstrng*(*)
      double precision x(nsd,nen),xs(nsd,nsd),det,shj(nsd+1,nen)
c
c...  local variables
c
cdebug      integer idb,jdb
c
cdebug      write(6,*) "Hello from getjac_f!"
c
      ierr=izero
c
c...calculate jacobian matrix for (x,y,z) to (r,s,t) transformation
c
cdebug      call dgemm("n","t",nsd,nsd,nen,one,shj,nsd+1,x,nsd,zero,xs,nsd)
cdebug      do idb=1,nen
cdebug        write(6,*) idb,(x(jdb,idb),jdb=1,nsd)
cdebug      end do
      call dgemm("n","t",nsd,nsd,nen,one,x,nsd,shj,nsd+ione,zero,xs,nsd)
c
c...form determinant of jacobian matrix and check for error condition
c
      det=xs(1,1)*xs(2,2)*xs(3,3)+xs(1,2)*xs(2,3)*xs(3,1)+xs(1,3)
     & *xs(2,1)*xs(3,2)-xs(1,3)*xs(2,2)*xs(3,1)-xs(1,2)*xs(2,1)
     & *xs(3,3)-xs(1,1)*xs(2,3)*xs(3,2)
cdebug      write(6,*) "iel,det:",iel,det
cdebug      do idb=1,nen
cdebug        write(6,*) "idb,x:",idb,(x(jdb,idb),jdb=1,nsd)
cdebug      end do
cdebug      do idb=1,nsd
cdebug        write(6,*) "idb,xs:",idb,(xs(jdb,idb),jdb=1,nsd)
cdebug      end do
cdebug      do idb=1,nen
cdebug        write(6,*) "idb,shj:",idb,(shj(jdb,idb),jdb=1,nsd+1)
cdebug      end do
cdebug      call flush(6)
      if(det.le.zero) then
        ierr=113
cdebug        write(6,700) iel,det
        write(errstrng,700) iel,det
      end if
c
 700  format("getjac:  element # ",i7,2x,1pe15.8)
      return
      end
c
c version
c $Id: getjac.f,v 1.8 2005/03/19 01:49:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
