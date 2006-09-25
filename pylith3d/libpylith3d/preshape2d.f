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
      subroutine preshape2d(sh,gauss,intord,ietype,nsnodes,nsgauss,ierr,
     & errstrng)
c
c...program to compute surface element integration information for a
c   given element type.
c
c   The arguments are:
c     Input:
c
c       intord                  = integration order
c                                 1 = full
c                                 2 = reduced
c                                 3 = Bbar
c       ietype                  = element type
c       nsnodes                 = number of nodes per element face
c       nsgauss                 = number of element face gauss points
c
c     Output:
c       sh(nsd,nsnodes,nsgauss)   = shape functions and their derivatives
c       gauss(nsd,nsgauss)        = Gauss point coordinates and weights
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
      integer intord,ietype,nsnodes,nsgauss,ierr
      double precision sh(nsd,nsnodes,nsgauss)
      double precision gauss(nsd,nsgauss)
      character*(*) errstrng
c
c...  local constants
c
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  user-defined functions
c
c
c...  local variables
c
      integer i,l,k,n,ind,nshsize,ngssize
      double precision rr,ss,drr,dss
      integer io(3)
c
cdebug      write(6,*) "Hello from preshape2d_f!"
c
c...  definitions
c
      ierr=izero
      ngssize=nsd*nsgauss
      nshsize=ngssize*nsnodes
c
c...  initialize arrays
c
      call fill(gauss,zero,ngssize)
      call fill(sh,zero,nshsize)
c
c... First type:  linear hex
c
      if(ietype.eq.1) then
        call plinhex2d(sh,gauss,nsnodes,nsgauss,intord)
c
c...  Types 2-27:  linear hex + infinite boundaries
c...  Type 28:  linear hex with one set of collapsed nodes (wrick)
c...  Type 29:  linear hex with two sets of collapsed nodes (wedge)
c...  Type 30:  linear hex with 4 points collapsed to a point (pyramid)
c     Leave these out for now.
c
      else if(ietype.lt.31) then
        ierr=106
        errstrng="preshape2d"
c
c...  Type 31:  linear tetrahedron
c     r, s, and t are tetrahedral coordinates.
c     One-point integration is used in all cases.
c
      else if(ietype.eq.31) then
        call plintet2d(sh,gauss,nsnodes,nsgauss,intord)
c
c... Type 32:  quadratic (20-node) hex
c...  Types 33-58:  quadratic hex + infinite boundaries.
c...  Type 59:  quad hex with one set of collapsed nodes (wrick).
c...  Type 60:  quad hex with three sets of collapsed nodes (wedge).
c...  Type 61:  quad hex with 9 points collapsed to a point (pyramid).
c...  Type 62  quadratic tetrahedron
c     Leave these out for now.
c
      else if(ietype.lt.62) then
        ierr=106
        errstrng="preshape"
      else
        ierr=106
        errstrng="preshape"
      end if
c
      return
      end
c
c version
c $Id: preshape.f,v 1.5 2005/03/22 02:22:18 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
