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
      subroutine assign_wink(winkdef,wink,iwinkdef,iwinkid,iwink,id,
     & numnp,nwink,nwinke)
c
c....program for reading and printing data on winkler restoring forces
c
c          winkdef(ndof,numnp) = values of winkler restoring spring
c                             constant, force(i,j)=-wink(i,j)*disp(i,j)
c
c          iwinkdef(ndof,numnp) = application mode:
c                               iwink = 0, no winkler forces,
c                               iwink = 1, applied throuthout computation
c                               iwink = -n, uses load history factor n
c
c          After assigning equation numbers, the following arrays are
c          returned:
c
c          wink(nwink)   = value of spring constant only for degrees of
c                          freedom with restoring forces.
c          wink(2,nwink) = application mode (1) and equation number (2)
c                          for each winkler restoring force.
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
      integer numnp,nwink,nwinke
      integer iwinkdef(ndof,nwinke),iwinkid(nwinke),iwink(2,nwink)
      integer id(ndof,numnp)
      double precision winkdef(ndof,nwinke),wink(nwink)
c
c...  included dimension and type statements
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n,nwtot,nnz
c
c...  included variable definitions
c
c
c...  loop over Winkler entries
c
      nwtot=izero
      do i=1,nwinke
        nnz=izero
        n=iwinkid(i)
        do j=1,ndof
          if(iwinkdef(j,i).ne.izero) then
            nnz=nnz+1
            nwtot=nwtot+1
            iwink(1,nwtot)=iwinkdef(j,i)
            iwink(2,nwtot)=id(j,n)
            wink(nwtot)=winkdef(j,i)
          end if
        end do
      end do
c
      return
      end
c
c version
c $Id: assign_wink.f,v 1.1 2005/04/16 00:35:38 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
