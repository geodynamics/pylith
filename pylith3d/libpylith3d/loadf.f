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
      subroutine loadf(fault,dfault,histry,deltp,nfault,nstep,
     & numfn,nhist,lastep,ierr,errstrng)
c
c...program to compute current split node displacements
c
c       if ihist = -1 constant velocity
c                =  0 constant displacement in elastic soln
c                =  # load history number
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
      integer nstep,numfn,nhist,lastep,ierr
      integer nfault(3,numfn)
      character errstrng*(*)
      double precision fault(ndof,numfn),dfault(ndof,numfn)
      double precision histry(nhist,lastep+1),deltp
c
c...  local variables
c
      integer l,j,ihist
cdebug      integer idb
c
cdebug      write(6,*) "Hello from loadf_f!"
c
      call fill(dfault,zero,ndof*numfn)
      do l=1,numfn
        ihist=nfault(3,l)
        if(ihist.gt.nhist) then
          ierr=100
          errstrng="loadf"
          return
        end if
        do j=1,ndof
          if(ihist.eq.-1) then
            dfault(j,l)=deltp*fault(j,l)
          else if((ihist.eq.izero).and.(nstep.eq.izero)) then
            dfault(j,l)=fault(j,l)
          else if(ihist.gt.izero) then
            dfault(j,l)=fault(j,l)*histry(ihist,nstep+1)
          end if
        end do
cdebug        write(6,"(3i7,3(2x,1pe15.8))") (nfault(idb,l),idb=1,3),
cdebug     &   (dfault(idb,l),idb=1,ndof)
      end do
      return
      end
c
c version
c $Id: loadf.f,v 1.3 2004/08/12 01:52:17 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
